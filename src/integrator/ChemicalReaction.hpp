/*
 Copyright (C) 2014
    Pierre de Buyl
    Jakub Krajniak (jkrajniak at gmail.com)
 Copyright (C) 2012,2013
    Max Planck Institute for Polymer Research
 Copyright (C) 2008,2009,2010,2011
    Max-Planck-Institute for Polymer Research & Fraunhofer SCAI

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

#include <utility>
#include <map>
#include <vector>

#include "types.hpp"
#include "logging.hpp"
#include "FixedPairList.hpp"
#include "Particle.hpp"
#include "SystemAccess.hpp"
#include "esutil/ESPPIterator.hpp"

#include "Extension.hpp"
#include "VelocityVerlet.hpp"
#include "VerletList.hpp"
#include "interaction/Potential.hpp"

#include "boost/signals2.hpp"

namespace espressopp {
namespace integrator {

/** Simple particle properties structure. Shortcut of Particle::ParticleProperties*/
typedef std::map<int, boost::shared_ptr<ParticleProperties> > TypeParticlePropertiesMap;

const int kCrCommTag = 0xad;

/** PostProcess **/

class PostProcess {
 public:
  PostProcess() { }
  virtual ~PostProcess() { }
  virtual bool operator()(Particle& p1, Particle& p2) = 0;

  /** Register this class so it can be used from Python. */
  static void registerPython();

 protected:
  shared_ptr<System> system_;
  static LOG4ESPP_DECL_LOGGER(theLogger);
};


class PostProcessChangesProperty : public integrator::PostProcess {
 public:
  bool operator()(Particle& p1, Particle& p2);
  void AddChangeProperty(int type_id, boost::shared_ptr<ParticleProperties> new_property);
  void RemoveChangeProperty(int type_id);

  /** Register this class so it can be used from Python. */
  static void registerPython();

 private:
  TypeParticlePropertiesMap type_properties_;
};


/** Abstract class for the chemical reactions. */
class Reaction {
 public:
  Reaction()
      : type_a_(-1),
        type_b_(-1),
        delta_a_(-1),
        delta_b_(-1),
        min_state_a_(-1),
        min_state_b_(-1),
        max_state_a_(-1),
        max_state_b_(-1),
        rate_(0.0),
        intramolecular_(false) {
    cutoff_ = 0.0;
    cutoff_sqr_ = 0.0;
  }

  Reaction(int type_a, int type_b, int delta_a, int delta_b, int min_state_a,
           int min_state_b, int max_state_a, int max_state_b, real cutoff,
           real rate, bool intramolecular = false)
      : type_a_(type_a),
        type_b_(type_b),
        delta_a_(delta_a),
        delta_b_(delta_b),
        min_state_a_(min_state_a),
        min_state_b_(min_state_b),
        max_state_a_(max_state_a),
        max_state_b_(max_state_b),
        rate_(rate),
        intramolecular_(intramolecular) {
    set_cutoff(cutoff);
  }
  virtual ~Reaction() { }

  void set_rate(real rate) { rate_ = rate; }
  real rate() { return rate_; }

  void set_cutoff(real cutoff) {
    cutoff_ = cutoff;
    cutoff_sqr_ = cutoff * cutoff;
  }
  real cutoff() { return cutoff_; }

  void set_type_a(int type_a) { type_a_ = type_a; }
  int type_a() { return type_a_; }

  void set_type_b(int type_b) { type_b_ = type_b; }
  int type_b() { return type_b_; }

  void set_delta_a(int delta_a) { delta_a_ = delta_a; }
  int delta_a() { return delta_a_; }

  void set_delta_b(int delta_b) { delta_b_ = delta_b; }
  int delta_b() { return delta_b_; }

  void set_min_state_a(int min_state_a) { min_state_a_ = min_state_a; }
  int min_state_a() { return min_state_a_; }

  void set_min_state_b(int min_state_b) { min_state_b_ = min_state_b; }
  int min_state_b() { return min_state_b_; }

  void set_max_state_a(int max_state_a) { max_state_a_ = max_state_a; }
  int max_state_a() { return max_state_a_; }

  void set_max_state_b(int max_state_b) { max_state_b_ = max_state_b; }
  int max_state_b() { return max_state_b_; }

  void set_intramolecular(bool intramolecular) { intramolecular_ = intramolecular; }
  bool intramolecular() { return intramolecular_; }

  void set_rng(const shared_ptr<esutil::RNG> rng) { rng_ = rng; }
  void set_interval(shared_ptr<int> interval) { interval_ = interval; }
  void set_dt(shared_ptr<real> dt) { dt_ = dt; }

  void AddPostProcess(const shared_ptr<integrator::PostProcess> pp) {
    post_process_.push_back(pp);
  }

  /** Checks if the pair is valid. */
  virtual bool IsValidPair(const Particle& p1, const Particle& p2);
  /** Checks if the pair has valid state. */
  virtual bool IsValidState(const Particle& p1, const Particle& p2);

  bool PostProcess(Particle &pA, Particle &pB);

  /** Register this class so it can be used from Python. */
  static void registerPython();

 protected:
  static LOG4ESPP_DECL_LOGGER(theLogger);

  int type_a_;  //!< type of reactant A
  int type_b_;  //!< type of reactant B
  int min_state_a_;  //!< min state of reactant A
  int min_state_b_;  //!< min state of reactant B
  int max_state_a_;  //!< max state of reactant A
  int max_state_b_;  //!< max state of reactant B
  int delta_a_;  //!< state change for reactant A
  int delta_b_;  //!< state change for reactant B
  real rate_;  //!< reaction rate
  real cutoff_;  //!< reaction cutoff
  real cutoff_sqr_;  //!< reaction cutoff^2

  bool intramolecular_;  //!< Allow to intramolecular reactions.

  shared_ptr<esutil::RNG> rng_;  //!< random number generator
  shared_ptr<int> interval_;  //!< number of steps between reaction loops
  shared_ptr<real> dt_;  //!< timestep from the integrator

  std::vector<shared_ptr<integrator::PostProcess> > post_process_;

};


/** Implements the synthesis reaction
 *
 * The reaction is in the form \f[ A + B \rightarrow A-B \f]
 *
 * In addtion the residue id of A molecule is transfer to B molecule so that
 * they have the same residue id. Other properties of A and B molecules remain the same.
 *
 */
class SynthesisReaction : public integrator::Reaction {
 public:
  SynthesisReaction(int type_a, int type_b, int delta_a, int delta_b,
                    int min_state_a, int min_state_b, int max_state_a,
                    int max_state_b, real cutoff, real rate,
                    bool intramolecular)
      : Reaction(type_a, type_b, delta_a, delta_b, min_state_a, min_state_b,
                 max_state_a, max_state_b, cutoff, rate, intramolecular) {
  }

  static void registerPython();
};

typedef boost::unordered_multimap<longint, std::pair<longint, int> > ReactionMap;
typedef std::vector<boost::shared_ptr<integrator::Reaction> > ReactionList;


/** Reaction scheme for polymer growth and curing/crosslinking

 This extension enables the rate-controlled stochastic curing of polymer
 systems, either for chain growth of step growth, depending on the
 parameters.

 The variables type_a, type_b, min_state_a, min_state_b, max_state_a, max_state_b
 control the particles that enter the curing reaction
 \f[ A^a + B^b \rightarrow A^{a+deltaA}-B^{b+deltaB} \f]
 where A and B may possess additional bonds not shown.

 An extra bond is added between A and B whenever the state of A and B falls
 into the defined range by variables min/max state.
 The condition is as follow:
 \f[ a >= minStateA \land stateA < maxStateA \f]
 the same holds for the particle B. Both condition should match.
 In addition if the intramolecular property is set to true (by default) then
 the reaction only happend between heterogenous molecules.

 The reaction proceeds by testing for all possible (A,B) pairs and
 selects them only at a given rate. It works in parallel, by gathering
 first the successful pairs between neigboring CPUs and ensuring that
 each particle enters only in one new bond per reaction step.
 */
class ChemicalReaction : public Extension {
 public:
  ChemicalReaction(shared_ptr<System> system,
                   shared_ptr<VerletList> _verletList,
                   shared_ptr<FixedPairList> _fixedPairList,
                   shared_ptr<storage::DomainDecomposition> _domdec);
  ~ChemicalReaction();

  void set_interval(int interval) {
    *interval_ = interval;
  }
  int interval() {
    return *interval_;
  }

  void Initialize();
  void AddReaction(boost::shared_ptr<integrator::Reaction> reaction);
  void RemoveReaction(int reaction_id);

  void React();

  void SendMultiMap(integrator::ReactionMap &mm);  //NOLINT
  void UniqueA(integrator::ReactionMap& potential_candidates);  //NOLINT
  void UniqueB(integrator::ReactionMap& potential_candidates,  //NOLINT
      integrator::ReactionMap& effective_candidates);  //NOLINT
  std::vector<Particle*> ApplyAR();

  /** Register this class so it can be used from Python. */
  static void registerPython();

 private:
  static LOG4ESPP_DECL_LOGGER(theLogger);
  void UpdateGhost(const std::vector<Particle*>& modified_particles);

  real current_cutoff_;

  shared_ptr<int> interval_;  //!< Number of steps between reaction loops.
  shared_ptr<real> dt_;  //!< Timestep from the integrator.

  shared_ptr<storage::DomainDecomposition> domdec_;
  shared_ptr<espressopp::interaction::Potential> potential_;
  shared_ptr<esutil::RNG> rng_;  //!< Random number generator.
  shared_ptr<FixedPairList> fixed_pair_list_;  //!< Bond list.
  shared_ptr<VerletList> verlet_list_;  //!< Verlet list of used potential

  boost::signals2::connection initialize_;
  boost::signals2::connection react_;

  integrator::ReactionMap potential_pairs_;  //!< Container for (A,B) potential partners/
  integrator::ReactionMap effective_pairs_;  //!< Container for (A,B) effective partners.

  ReactionList reaction_list_;  //<! Container for reactions.

  void connect();
  void disconnect();
};

}
  // namespace integrator
} // namespace espressopp

#endif
