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

/** Simple particle properties structure. Shortcut of Particle::ParticleProperties*/

typedef std::pair<Particle*, Particle*> ParticlePair;

const int kCrCommTag = 0xad;

/** Class for the chemical reactions. */
class Reaction {
 public:
  Reaction()
      : type_1_(-1),
        type_2_(-1),
        delta_1_(-1),
        delta_2_(-1),
        min_state_1_(-1),
        max_state_1_(-1),
        min_state_2_(-1),
        max_state_2_(-1),
        rate_(0.0),
        reverse_(false),
        intramolecular_(false), active_(true) {
    cutoff_ = 0.0;
    cutoff_sqr_ = 0.0;
    min_cutoff_ = 0.0;
    min_cutoff_sqr_ = 0.0;
  }

  Reaction(int type_1, int type_2, int delta_1, int delta_2, int min_state_1,
           int max_state_1, int min_state_2, int max_state_2, real cutoff,
           real rate, shared_ptr<FixedPairList> fpl,
           bool intramolecular = false)
      : type_1_(type_1),
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
        active_(true) {
    set_cutoff(cutoff);
    set_min_cutoff(0.0);
  }
  virtual ~Reaction() { }

  void set_rate(real rate) { rate_ = rate; }
  real rate() { return rate_; }

  void set_cutoff(real cutoff) {
    cutoff_ = cutoff;
    cutoff_sqr_ = cutoff * cutoff;
  }
  real cutoff() { return cutoff_; }

  void set_min_cutoff(real cutoff) {
    min_cutoff_ = cutoff;
    min_cutoff_sqr_ = cutoff * cutoff;
  }
  real min_cutoff() { return min_cutoff_; }

  void set_type_1(int type_1) { type_1_ = type_1; }
  int type_1() { return type_1_; }

  void set_type_2(int type_2) { type_2_ = type_2; }
  int type_2() { return type_2_; }

  void set_delta_1(int delta_1) { delta_1_ = delta_1; }
  int delta_1() { return delta_1_; }

  void set_delta_2(int delta_2) { delta_2_ = delta_2; }
  int delta_2() { return delta_2_; }

  void set_min_state_1(int min_state_1) { min_state_1_ = min_state_1; }
  int min_state_1() { return min_state_1_; }

  void set_min_state_2(int min_state_2) { min_state_2_ = min_state_2; }
  int min_state_2() { return min_state_2_; }

  void set_max_state_1(int max_state_1) { max_state_1_ = max_state_1; }
  int max_state_1() { return max_state_1_; }

  void set_max_state_2(int max_state_2) { max_state_2_ = max_state_2; }
  int max_state_2() { return max_state_2_; }

  void set_intramolecular(bool intramolecular) { intramolecular_ = intramolecular; }
  bool intramolecular() { return intramolecular_; }

  void set_rng(const shared_ptr<esutil::RNG> rng) { rng_ = rng; }
  void set_interval(shared_ptr<int> interval) { interval_ = interval; }
  void set_dt(shared_ptr<real> dt) { dt_ = dt; }
  void set_bc(bc::BC *bc) { bc_ = bc; }

  bool reverse() { return reverse_; };
  void set_reverse(bool r) { reverse_ = r; }

  bool active() { return active_; }
  void set_active(bool s) { active_ = s; }

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

  /** Checks if the pair is valid. */
  virtual bool IsValidPair(Particle& p1, Particle& p2, ParticlePair &correct_order);
  /** Checks if the pair has valid state. */
  bool IsValidState(Particle& p1, Particle& p2, ParticlePair &correct_order);

  bool IsValidStateT_1(Particle &p);
  bool IsValidStateT_2(Particle &p);

  std::set<Particle*> PostProcess_T1(Particle &p, Particle &partner);
  std::set<Particle*> PostProcess_T2(Particle &p, Particle &partner);

  shared_ptr<FixedPairList> fixed_pair_list_;  //!< Bond list.

  /** Register this class so it can be used from Python. */
  static void registerPython();

 protected:
  static LOG4ESPP_DECL_LOGGER(theLogger);

  int type_1_;  //!< type of reactant A
  int type_2_;  //!< type of reactant B
  int min_state_1_;  //!< min state of reactant A
  int min_state_2_;  //!< min state of reactant B
  int max_state_1_;  //!< max state of reactant A
  int max_state_2_;  //!< max state of reactant B
  int delta_1_;  //!< state change for reactant A
  int delta_2_;  //!< state change for reactant B
  real rate_;  //!< reaction rate
  real cutoff_;  //!< reaction cutoff
  real cutoff_sqr_;  //!< reaction cutoff^2
  real min_cutoff_;  //!< min reaction cutoff
  real min_cutoff_sqr_;  //!< min reaction cutoff^2
  bool active_ ;  //!< is reaction active, by default true

  bool intramolecular_;  //!< Allow to intramolecular reactions.

  bool reverse_;  //!< If true then reaction will break a bond.

  shared_ptr<esutil::RNG> rng_;  //!< random number generator
  shared_ptr<int> interval_;  //!< number of steps between reaction loops
  shared_ptr<real> dt_;  //!< timestep from the integrator
  bc::BC *bc_;  //!< boundary condition

  std::vector<shared_ptr<integrator::ChemicalReactionPostProcess> > post_process_T1;
  std::vector<shared_ptr<integrator::ChemicalReactionPostProcess> > post_process_T2;
};


class DissociationReaction : public Reaction {
 public:
  DissociationReaction() : Reaction() {
    cutoff_ = 0.0;
    cutoff_sqr_ = 0.0;
    reverse_ = true;
  }

  DissociationReaction(
      int type_1, int type_2, int delta_1, int delta_2, int min_state_1, int max_state_1,
      int min_state_2, int max_state_2,
      real break_cutoff,
      real break_rate,
      shared_ptr<FixedPairList> fpl
      ) : Reaction(type_1, type_2, delta_1, delta_2,
          min_state_1, max_state_1, min_state_2, max_state_2, break_cutoff, break_rate, fpl,
          true), diss_rate_(0.0) {
    reverse_ = true;
  }
  virtual ~DissociationReaction() { }

  real diss_rate() { return diss_rate_; }
  void set_diss_rate(real s) { diss_rate_ = s; }

  bool IsValidPair(Particle& p1, Particle& p2, ParticlePair &correct_order);

  /** Register this class so it can be used from Python. */
  static void registerPython();

 protected:
  static LOG4ESPP_DECL_LOGGER(theLogger);

 private:
  real diss_rate_;  //!< Dissociation rate.

};

}  // namespace integrator
}  // namespace espressopp

#endif
