/*
  Copyright (C) 2014
      Pierre de Buyl, Jakub Krajniak
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


namespace espresso {
namespace integrator {

const int kCrCommTag = 0xad;

/** Reaction scheme for polymer growth and curing/crosslinking

This extension enables the rate-controlled stochastic curing of polymer
systems, either for chain growth of step growth, depending on the
parameters.

The variables type_a, type_b, min_state_a, min_state_b, max_state_a, max_state_b
control the particles that enter the curing reaction
\f[ A^a + B^0 -> A^{a+deltaA}-B^{deltaB} \f]
where A and B may possess additional bonds not shown.

An extra bond is added between A and B whenever the state of A and B falls
into the defined range by variables min/max state.
The condition is as follow: \f[state_A >= min_state_a & state_A <= max_state_a \f]
the same holds for the particle B.

The reaction proceeds by testing for all possible (A,B) pairs and
selects them only at a given rate. It works in parallel, by gathering
first the successful pairs between neigboring CPUs and ensuring that
each particle enters only in one new bond per reaction step.
*/

class Reaction {
 public:
  Reaction();
  virtual ~Reaction() {}

  void set_rate(real rate) { rate_ = rate; }
  real rate() { return rate_; }

  void set_cutoff(real cutoff) {
      cutoff_ = cutoff;
      cutoff_sqr_ = cutoff * cutoff;
  }
  real cutoff() { return cutoff_; }

  void set_type_a(int typeA) { type_a_ = typeA; }
  int type_a() { return type_a_; }

  void set_type_b(int typeB) { type_b_ = typeB; }
  int type_b() { return type_b_; }

  void set_delta_a(int deltaA) { delta_a_ = deltaA; }
  int delta_a() { return delta_a_; }

  void set_delta_b(int deltaB) { delta_b_ = deltaB; }
  int delta_b() { return delta_b_; }

  void set_min_state_a(int minStateA) { min_state_a_ = minStateA; }
  int min_state_a() { return min_state_a_; }

  void set_min_state_b(int minStateB) { min_state_b_ = minStateB; }
  int min_state_b() { return min_state_b_; }

  void set_max_state_a(int maxStateA) { max_state_a_ = maxStateA; }
  int max_state_a() { return max_state_a_; }

  void set_max_state_b(int maxStateB) { max_state_b_ = maxStateB; }
  int max_state_b() { return max_state_b_; }

  void set_rng(const shared_ptr<esutil::RNG> rng) { rng_ = rng; }
  void set_interval(shared_ptr<int> interval) { interval_ = interval; }
  void set_dt(shared_ptr<real> dt) { dt_ = dt; }

  virtual bool IsValidPair(const Particle& p1, const Particle& p2);

  /** Register this class so it can be used from Python. */
  static void registerPython();

 protected:
  int type_a_;  //!< type of reactant A
  int type_b_;  //!< type of reactant B
  int min_state_a_;  //!< min state of reactant A
  int min_state_b_;  //!< min state of reactant B
  int max_state_a_;  //!< max state of reactant A
  int max_state_b_;  //!< max state of reactant B
  int delta_a_;  //!< state change for reactant A
  int delta_b_;  //!< state change for reactant B
  int rate_;  //!< reaction rate
  real cutoff_;  //!< reaction cutoff
  real cutoff_sqr_;  //!< reactio cutoff^2

  shared_ptr<esutil::RNG> rng_;  //!< random number generator
  shared_ptr<int> interval_;  //!< number of steps between reaction loops
  shared_ptr<real> dt_;  //!< timestep from the integrator
};

class SynthesisReaction : public Reaction {
 public:
  bool IsValidPair(const Particle& p1, const Particle &p2);
  static void registerPython();
};


typedef boost::unordered_multimap<longint, std::pair<longint, longint> > ReactionMap;
typedef std::vector<shared_ptr<integrator::Reaction> > ReactionList;


class ChemicalReaction:public Extension {
 public:
  ChemicalReaction(shared_ptr<System> system,
                   shared_ptr<VerletList> _verletList,
                   shared_ptr<FixedPairList> _fixedPairList,
                   shared_ptr<storage::DomainDecomposition> _domdec);
  ~ChemicalReaction();

  void set_interval(int interval) { *interval_ = interval; }
  int interval() { return *interval_; }

  void Initialize();
  void AddReaction(shared_ptr<integrator::Reaction> reaction);
  void RemoveReaction(int reaction_id);

  void React();

  void SendMultiMap(integrator::ReactionMap &mm);  //NOLINT
  void UniqueA(integrator::ReactionMap& mm);  //NOLINT
  void UniqueB(integrator::ReactionMap& mm,  //NOLINT
               integrator::ReactionMap& nn);  //NOLINT
  void ApplyAR();

  /** Register this class so it can be used from Python. */
  static void registerPython();

 private:
  static LOG4ESPP_DECL_LOGGER(theLogger);

  real current_cutoff_;

  shared_ptr<int> interval_;  //!< number of steps between reaction loops
  shared_ptr<real> dt_;  //!< timestep from the integrator

  shared_ptr<storage::DomainDecomposition> domdec_;
  shared_ptr<espresso::interaction::Potential> potential_;
  shared_ptr<esutil::RNG> rng_;  //!< random number generator
  shared_ptr<FixedPairList> fixed_pair_list_;
  shared_ptr<VerletList> verlet_list_;

  boost::signals2::connection initialize_;
  boost::signals2::connection react_;

  /** container for (A,B) potential partners */
  integrator::ReactionMap potential_pairs_;
  /** container for (A,B) effective partners */
  integrator::ReactionMap effective_pairs_;

  /** container for reactions */
  ReactionList reaction_list_;

  void connect();
  void disconnect();
};

}  // namespace integrator
}  // namespace espresso

#endif
