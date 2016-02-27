/*
   Copyright (C) 2014-2016
   Jakub Krajniak (jkrajniak at gmail.com)

   This file is part of ESPResSo++.

   ESPResSo++ is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   ESPResSo++ is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

// ESPP_CLASS
#ifndef _INTEGRATOR_CHEMICALREACTION_HPP
#define _INTEGRATOR_CHEMICALREACTION_HPP

#include <boost/unordered_map.hpp>
#include <boost/signals2.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/variate_generator.hpp>

#include <utility>
#include <map>
#include <set>
#include <vector>

#include "types.hpp"
#include "logging.hpp"
#include "FixedPairList.hpp"
#include "Particle.hpp"
#include "SystemAccess.hpp"
#include "esutil/ESPPIterator.hpp"

#include "integrator/Extension.hpp"
#include "integrator/VelocityVerlet.hpp"
#include "VerletList.hpp"
#include "interaction/Potential.hpp"

#include "ChemicalReactionPostProcess.hpp"

namespace espressopp {
namespace integrator {
typedef std::pair<Particle *, Particle *> ParticlePair;

const int kCrCommTag = 0xad;// @warning: this made problems when multiple extension will be enabled.

class ReactionCutoff {
public:
  ReactionCutoff() { }

  virtual ~ReactionCutoff() { }

  virtual bool check(Particle &p1, Particle &p2) = 0;
  virtual real cutoff() = 0;

  /** Register this class so it can be used from Python. */
  static void registerPython();

protected:
  static LOG4ESPP_DECL_LOGGER(theLogger);
};

class ReactionCutoffStatic : public ReactionCutoff {
public:
  ReactionCutoffStatic() {
    set_cutoff(0.0, 0.0);
  }

  ReactionCutoffStatic(real min_cutoff, real max_cutoff) {
    set_cutoff(min_cutoff, max_cutoff);
  }

  bool check(Particle &p1, Particle &p2);
  real cutoff() {return max_cutoff_; }

  /** Register this class so it can be used from Python. */
  static void registerPython();

protected:
  static LOG4ESPP_DECL_LOGGER(theLogger);

private:
  void set_cutoff(real min_cutoff, real max_cutoff) {
    min_cutoff_ = min_cutoff;
    min_cutoff_sqr_ = min_cutoff * min_cutoff;
    max_cutoff_ = max_cutoff;
    max_cutoff_sqr_ = max_cutoff * max_cutoff;
  }

  real min_cutoff_;
  real min_cutoff_sqr_;
  real max_cutoff_;
  real max_cutoff_sqr_;
};

class ReactionCutoffRandom : public ReactionCutoff {
public:
  ReactionCutoffRandom(real eq_distance, real eq_width, longint seed)
        :eq_width_(eq_width), eq_distance_(eq_distance), seed_(seed),
            generator_(boost::mt19937(seed), boost::normal_distribution<>(0.0, eq_width)) { }

  bool check(Particle &p1, Particle &p2);
  real cutoff() {return eq_distance_ + 0.5 * eq_width_; }

  /** Register this class so it can be used from Python. */
  static void registerPython();

protected:
  static LOG4ESPP_DECL_LOGGER(theLogger);

private:

  real eq_distance_;
  real eq_width_;
  longint interval_;
  longint seed_;
  boost::variate_generator<boost::mt19937, boost::normal_distribution<> > generator_;
};

/** Base class for performing addition reaction.
 *
 * It modeled following reaction:
 *
 *  \f[ A^a + B^b \rightarrow A^{a+deltaA}-B^{b+deltaB} \f]
 *
 **/
class Reaction {
public:
  Reaction()
        :type_1_(-1),
            type_2_(-1),
            delta_1_(-1),
            delta_2_(-1),
            min_state_1_(-1),
            max_state_1_(-1),
            min_state_2_(-1),
            max_state_2_(-1),
            rate_(0.0),
            reverse_(false),
            intramolecular_(false), active_(true) { }

  /*** Constructor of Reaction object.
   *
   * \f[ A^a + B^b \rightarrow A^{a+deltaA}-B^{b+deltaB} \f]
   *
   * @param type_1 The type of A particle.
   * @param type_2 The type of B particle.
   * @param delta_1 The delta A.
   * @param delta_2 The delta B.
   * @param min_state_1 The minimum state of particle A (greater equal to it).
   * @param max_state_1 The maximum state of particle A (less not equal to it).
   * @param min_state_2 The minimum state of particle B.
   * @param max_state_2 The maximum state of particle B.
   * @param cutoff The reaction cutoff distance.
   * @param rate The reaction rate.
   * @param fpl The espressopp.FixedPairList with the new bonds that are added by reaction.
   * @param intramolecular If set to true then intramolecular bonds are allowed.
   *
   */
  Reaction(int type_1, int type_2, int delta_1, int delta_2, int min_state_1, int max_state_1, int
      min_state_2, int max_state_2, real rate, shared_ptr<FixedPairList> fpl, bool
      intramolecular = false)
        :type_1_(type_1),
            type_2_(type_2),
            delta_1_(delta_1),
            delta_2_(delta_2),
            min_state_1_(min_state_1),
            max_state_1_(max_state_1),
            min_state_2_(min_state_2),
            max_state_2_(max_state_2),
            rate_(rate),
            reverse_(false),
            fixed_pair_list_(fpl),
            intramolecular_(intramolecular),
            active_(true) { }

  virtual ~Reaction() { }

  virtual real cutoff() {
    return reaction_cutoff_->cutoff();
  }

  void set_rate(real rate) {rate_ = rate; }

  real rate() {return rate_; }

  void set_type_1(int type_1) {type_1_ = type_1; }

  int type_1() {return type_1_; }

  void set_type_2(int type_2) {type_2_ = type_2; }

  int type_2() {return type_2_; }

  void set_delta_1(int delta_1) {delta_1_ = delta_1; }

  int delta_1() {return delta_1_; }

  void set_delta_2(int delta_2) {delta_2_ = delta_2; }

  int delta_2() {return delta_2_; }

  void set_min_state_1(int min_state_1) {min_state_1_ = min_state_1; }

  int min_state_1() {return min_state_1_; }

  void set_min_state_2(int min_state_2) {min_state_2_ = min_state_2; }

  int min_state_2() {return min_state_2_; }

  void set_max_state_1(int max_state_1) {max_state_1_ = max_state_1; }

  int max_state_1() {return max_state_1_; }

  void set_max_state_2(int max_state_2) {max_state_2_ = max_state_2; }

  int max_state_2() {return max_state_2_; }

  void set_intramolecular(bool intramolecular) {intramolecular_ = intramolecular; }

  bool intramolecular() {return intramolecular_; }

  void set_rng(const shared_ptr<esutil::RNG> rng) {rng_ = rng; }

  void set_interval(shared_ptr<int> interval) {interval_ = interval; }

  void set_dt(shared_ptr<real> dt) {dt_ = dt; }

  void set_bc(bc::BC *bc) {bc_ = bc; }

  bool reverse() {return reverse_; };

  void set_reverse(bool r) {reverse_ = r; }

  bool active() {return active_; }

  /** Activate the reaction*/
  void set_active(bool s) {active_ = s; }

  /**
   * Define post-process method after reaction occures.
   *
   * @param pp: The PostProcess object.
   * @param type: if set to 0 then applied to both particle types,
   *     if set to 1 then applied only to type_1 particle,
   *     if set to 2 then applied only to type_2 particle.
   */
  void AddPostProcess(const shared_ptr<integrator::ChemicalReactionPostProcess> pp, int type = 0) {
    switch (type) {
      case 1:
        post_process_T1.push_back(pp); break;
      case 2:
        post_process_T2.push_back(pp); break;
      case 0:
        post_process_T1.push_back(pp);
        post_process_T2.push_back(pp);
        break;
    }
  }

  void SetReactionCutoff(shared_ptr<ReactionCutoff> rc) {
    reaction_cutoff_ = rc;
  }

  /** Checks if the pair is valid. */
  virtual bool IsValidPair(Particle &p1, Particle &p2, ParticlePair &correct_order);

  /** Checks if the pair has valid state. */
  bool IsValidState(Particle &p1, Particle &p2, ParticlePair &correct_order);

  bool IsValidStateT_1(Particle &p);
  bool IsValidStateT_2(Particle &p);

  std::set<Particle *> PostProcess_T1(Particle &p, Particle &partner);
  std::set<Particle *> PostProcess_T2(Particle &p, Particle &partner);

  shared_ptr<FixedPairList> fixed_pair_list_;//!< Bond list.

  /** Register this class so it can be used from Python. */
  static void registerPython();

protected:
  static LOG4ESPP_DECL_LOGGER(theLogger);

  int type_1_;//!< type of reactant A
  int type_2_;//!< type of reactant B
  int min_state_1_;//!< min state of reactant A
  int min_state_2_;//!< min state of reactant B
  int max_state_1_;//!< max state of reactant A
  int max_state_2_;//!< max state of reactant B
  int delta_1_;//!< state change for reactant A
  int delta_2_;//!< state change for reactant B
  real rate_;//!< reaction rate
  bool active_;//!< is reaction active, by default true

  bool intramolecular_;//!< Allow to intramolecular reactions.

  bool reverse_;//!< If true then reaction will break a bond.

  shared_ptr<esutil::RNG> rng_;//!< random number generator
  shared_ptr<int> interval_;//!< number of steps between reaction loops
  shared_ptr<real> dt_;//!< timestep from the integrator

  bc::BC *bc_;//!< boundary condition

  std::vector<shared_ptr<integrator::ChemicalReactionPostProcess> > post_process_T1;
  std::vector<shared_ptr<integrator::ChemicalReactionPostProcess> > post_process_T2;

  shared_ptr<ReactionCutoff> reaction_cutoff_;
};

/*** Defines dissociation reactions.
 *
 * This special type of reaction, modeling follwoing reaction:
 *
 * \f[ A:B -> A + B \f]
 *
 * When this reaction is invoked, the list of bonds between particles
 * A and B is scanned, the bond is removed when those conditions occures:
 *  - distance between A and B is larger than specific cut_off distance,
 *  - the \f[ k\Delta t \Phi < W\f] when \f[k\f] is a rate
 *
 * It is also possible to remove the bond without checking if the distance exceed the bond
 * by defining only the diss_rate. By default, this rate is set to 0
 */
class DissociationReaction : public Reaction {
public:
  DissociationReaction():
    Reaction() {
    break_cutoff_ = 0.0;
    break_cutoff_sqr_ = 0.0;
    reverse_ = true;
  }

  DissociationReaction(int type_1, int type_2, int delta_1, int delta_2, int min_state_1, int
      max_state_1, int min_state_2, int max_state_2, real break_cutoff, real break_rate,
      shared_ptr<FixedPairList> fpl):
    Reaction(type_1, type_2, delta_1, delta_2,
        min_state_1, max_state_1, min_state_2, max_state_2, break_rate, fpl,
        true), diss_rate_(0.0), break_cutoff_(break_cutoff) {
    reverse_ = true;
    break_cutoff_sqr_ = break_cutoff_ * break_cutoff_;
  }

  virtual ~DissociationReaction() { }

  real diss_rate() {return diss_rate_; }

  void set_diss_rate(real s) {diss_rate_ = s; }

  void set_cutoff(real cutoff) {break_cutoff_ = cutoff; break_cutoff_sqr_ = cutoff * cutoff; }

  real cutoff() {return break_cutoff_; }

  bool IsValidPair(Particle &p1, Particle &p2, ParticlePair &correct_order);

  /** Register this class so it can be used from Python. */
  static void registerPython();

protected:
  static LOG4ESPP_DECL_LOGGER(theLogger);

private:

  real diss_rate_;//!< Dissociation rate.
  real break_cutoff_;
  real break_cutoff_sqr_;
};
}// namespace integrator
}// namespace espressopp

#endif
