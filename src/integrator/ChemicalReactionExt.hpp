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
   along with this program.  If not, see <http:  //www.gnu.org/licenses/>.
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
#include "esutil/Timer.hpp"

#include "integrator/Extension.hpp"
#include "integrator/VelocityVerlet.hpp"
#include "VerletList.hpp"
#include "interaction/Potential.hpp"

#include "ChemicalReaction.hpp"
#include "ChemicalReactionPostProcess.hpp"

namespace espressopp {
namespace integrator {

struct ReactionDef {
  longint reaction_id;
  real reaction_rate;
  real reaction_r_sqr;
  // 1 -> first => T1, second.first => T2
  // 2 -> first => T2, second.first => T1
  uint8_t order;

  ReactionDef(longint r_id, real rr, real r_sqr) {
    reaction_id = r_id;
    reaction_rate = rr;
    reaction_r_sqr = r_sqr;
  }
};

typedef std::multimap<longint, std::pair<longint, ReactionDef> > ReactionMap;
typedef std::vector<boost::shared_ptr<integrator::Reaction> > ReactionList;

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
   The conditions are as follow:
   \f[ a >= minStateA \land stateA < maxStateA \f]
   the same holds for the particle B. Both condition should match.

   In addition to that, there is also a probabilistic condition
   \f[ k\Delta t \Phi < W \f]
   where \f$k\f$ is a kinetic rate, \f$\Delta t\f$ is an integrator time step,
   \f$W\f$ is a number from uniform random number generator and \f$\Phi\f$ is
   an interval between invokation of the reactions.

   It is possible to exclude reactions between particles that are in the same
   molecule (like polymer chain etc.). This can be done by settings intramolecular
   property to false (default).

   The reaction proceeds by testing for all possible (A,B) pairs and
   selects them only at a given rate. It works in parallel, by gathering
   first the successful pairs between neigboring CPUs and ensuring that
   each particle enters only in one new bond per reaction step.
 */

class ChemicalReaction : public Extension {
public:
  ChemicalReaction(shared_ptr<System> system, shared_ptr<VerletList> _verletList,
      shared_ptr<storage::DomainDecomposition> _domdec, shared_ptr<TopologyManager> tm);
  ~ChemicalReaction();

  void set_interval(int interval) {
    *interval_ = interval;
  }

  /// Gets interval between when reaction happens.
  int interval() {
    return *interval_;
  }

  bool is_nearest() { return is_nearest_; }
  void set_is_nearest(bool s_) { is_nearest_ = s_; }

  longint bond_limit() { return bond_limit_; }
  void set_bond_limit(longint s_) { bond_limit_ = s_; }

  /** Register this class so it can be used from Python. */
  static void registerPython();

private:
  static LOG4ESPP_DECL_LOGGER(theLogger);
  void addReaction(boost::shared_ptr<integrator::Reaction> reaction);

  void React();

  void sendMultiMap(integrator::ReactionMap &mm);
  void UniqueA(integrator::ReactionMap &potential_candidates);
  void UniqueB(integrator::ReactionMap &potential_candidates,
      integrator::ReactionMap &effective_candidates);
  void ApplyAR(std::set<Particle *> &modified_particles);
  void ApplyDR(std::set<Particle *> &modified_particles);

  void updateGhost(const std::set<Particle *> &modified_particles);

  real current_cutoff_;  //!< Maximal cutoff use for VerletList.

  shared_ptr<int> interval_;  //!< Number of steps between reaction loops.
  shared_ptr<real> dt_;  //!< Timestep from the integrator.

  shared_ptr<storage::DomainDecomposition> domdec_;
  shared_ptr<espressopp::interaction::Potential> potential_;
  shared_ptr<esutil::RNG> rng_;  //!< Random number generator.
  shared_ptr<VerletList> verlet_list_;  //!< Verlet list of used potential

  boost::signals2::connection react_;

  integrator::ReactionMap potential_pairs_;  //!< Container for (A,B) potential partners/
  integrator::ReactionMap effective_pairs_;  //!< Container for (A,B) effective partners.

  ReactionList reaction_list_;  // <! Container for reactions.
  ReactionList reverse_reaction_list_;  // <! Container for reverse reactions.

  shared_ptr<TopologyManager> tm_;  //<! TopologyManager object.

  bool is_nearest_;  //!< If set to True then nearest neighbour is taken instead of random particle.

  longint bond_limit_;

  void connect();
  void disconnect();

  ///Timers
  esutil::WallTimer wallTimer;  //!< used for timing

  real timeComm;
  real timeUpdateGhost;
  real timeApplyAR;
  real timeApplyDR;
  real timeLoopPair;

  void resetTimers() {
    timeComm = 0.0;
    timeUpdateGhost = 0.0;
    timeApplyAR = 0.0;
    timeApplyDR = 0.0;
    timeLoopPair = 0.0;
  }

  python::list getTimers();

  void printMultiMap(ReactionMap &rmap, std::string comment);
};
}  // namespace integrator
}  // namespace espressopp
#endif
