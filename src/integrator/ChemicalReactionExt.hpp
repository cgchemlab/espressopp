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
#ifndef _INTEGRATOR_CHEMICALREACTIONEXT_HPP
#define _INTEGRATOR_CHEMICALREACTIONEXT_HPP

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

#include "ChemicalReaction.hpp"
#include "ChemicalReactionPostProcess.hpp"


namespace espressopp {
namespace integrator {


typedef boost::unordered_multimap<longint, std::pair<longint, int> > ReactionMap;
typedef std::vector<boost::shared_ptr<integrator::Reaction> > ReactionList;
typedef std::vector<boost::unordered_multimap<longint, longint> > RevReactionPairList;


/** Reaction scheme for polymer growth and curing/crosslinking

 This extension enables the rate-controlled stochastic curing of polymer
 systems, either for chain growth of step growth, depending on the
 parameters.

 The variables type_1, type_2, min_state_1, min_state_2, max_state_1, max_state_2
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

  void React();

  void SendMultiMap(integrator::ReactionMap &mm);  //NOLINT
  void UniqueA(integrator::ReactionMap& potential_candidates);  //NOLINT
  void UniqueB(integrator::ReactionMap& potential_candidates,  //NOLINT
      integrator::ReactionMap& effective_candidates);  //NOLINT
  void ApplyAR(std::set<Particle*>& modified_particles);
  void ApplyDR(std::set<Particle*>& modified_particles);

  /** Register this class so it can be used from Python. */
  static void registerPython();

 private:
  static LOG4ESPP_DECL_LOGGER(theLogger);
  void UpdateGhost(const std::set<Particle*>& modified_particles);

  real current_cutoff_;

  shared_ptr<int> interval_;  //!< Number of steps between reaction loops.
  shared_ptr<real> dt_;  //!< Timestep from the integrator.

  shared_ptr<storage::DomainDecomposition> domdec_;
  shared_ptr<espressopp::interaction::Potential> potential_;
  shared_ptr<esutil::RNG> rng_;  //!< Random number generator.
  shared_ptr<VerletList> verlet_list_;  //!< Verlet list of used potential

  boost::signals2::connection initialize_;
  boost::signals2::connection react_;

  integrator::ReactionMap potential_pairs_;  //!< Container for (A,B) potential partners/
  integrator::ReactionMap effective_pairs_;  //!< Container for (A,B) effective partners.

  ReactionList reaction_list_;  //<! Container for reactions.
  ReactionList reverse_reaction_list_;  //<! Container for reverse reactions.

  void connect();
  void disconnect();
};


}  // namespace integrator
}  // namespace espressopp

#endif
