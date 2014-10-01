/*
  Copyright (C) 2014
      Pierre de Buyl
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

#include "ChemicalReaction.hpp"

#include <utility>
#include <vector>
#include "python.hpp"

#include "types.hpp"
#include "System.hpp"
#include "storage/Storage.hpp"
#include "iterator/CellListIterator.hpp"
#include "esutil/RNG.hpp"
#include "bc/BC.hpp"
#include "storage/NodeGrid.hpp"
#include "storage/DomainDecomposition.hpp"

#include "boost/make_shared.hpp"


namespace espresso {
namespace integrator {

LOG4ESPP_LOGGER(ChemicalReaction::theLogger, "ChemicalReaction");

Reaction::Reaction():rate_(0.0),
                     type_a_(-1),
                     type_b_(-1),
                     delta_a_(-1),
                     delta_b_(-1),
                     min_state_a_(0),
                     min_state_b_(0),
                     max_state_a_(0),
                     max_state_b_(0),
                     cutoff_(0.0),
                     cutoff_sqr_(0.0) { }


bool SynthesisReaction::IsValidPair(const Particle& p1, const Particle& p2) {
  Real3D distance = p1.position() - p2.position();
  real distance_2 = distance.sqr();

  if ((distance_2 < cutoff_sqr_) && ((*rng_)() < rate_*(*dt_)*(*interval_))) {
    int p1_state = p1.state();
    int p2_state = p2.state();
    if ((p1.type() == type_a_) && (p2.type() == type_b_) &&
        (p1_state >= min_state_a_ && p1_state <= max_state_a_) &&
        (p2_state >= min_state_b_ && p2_state <= max_state_b_)) {
      return true;
    } else if ((p1.type() == type_b_) && (p2.type() == type_a_) &&
        (p1_state >= min_state_b_ && p1_state <= max_state_b_) &&
        (p2_state >= min_state_a_ && p2_state <= max_state_a_)) {
      return true;
    }
  }
  return false;
}

void SynthesisReaction::registerPython() {
  using namespace espresso::python;  //NOLINT
  class_<SynthesisReaction, shared_ptr<SynthesisReaction> >
    ("integrator_SynthesisReaction", no_init)
      .add_property(
        "type_a",
        &Reaction::type_a,
        &Reaction::set_type_a)
      .add_property(
        "type_b",
        &Reaction::type_b,
        &Reaction::set_type_b)
      .add_property(
        "delta_a",
        &Reaction::delta_a,
        &Reaction::set_delta_a)
      .add_property(
        "min_state_a",
        &Reaction::min_state_a,
        &Reaction::set_min_state_a)
      .add_property(
        "max_state_a",
        &Reaction::max_state_b,
        &Reaction::set_max_state_b)
      .add_property(
        "delta_b",
        &Reaction::delta_b,
        &Reaction::set_delta_b)
      .add_property(
        "min_state_b",
        &Reaction::min_state_b,
        &Reaction::set_min_state_b)
      .add_property(
        "max_state_b",
        &Reaction::max_state_b,
        &Reaction::set_max_state_b)
      .add_property(
        "rate",
        &Reaction::rate,
        &Reaction::set_rate)
      .add_property(
        "cutoff",
        &Reaction::cutoff,
        &Reaction::set_cutoff);
  }


/** ChemicalReaction */
ChemicalReaction::ChemicalReaction(shared_ptr<System> system,
                                   shared_ptr<VerletList> verletList,
                                   shared_ptr<FixedPairList> fpl,
                                   shared_ptr<storage::DomainDecomposition> domdec)
                 :
                 Extension(system),
                 verlet_list_(verletList),
                 fixed_pair_list_(fpl),
                 domdec_(domdec) {
  type = Extension::Reaction;

  current_cutoff_ = verletList->getVerletCutoff() - system->getSkin();

  if (!system->rng) throw std::runtime_error("System has no RNG.");

  rng_ = system->rng;
  LOG4ESPP_INFO(theLogger, "ChemicalReaction constructed");
  dt_ = boost::make_shared<real>();
  interval_ = boost::make_shared<int>();

  reaction_list_ = ReactionList();
}


ChemicalReaction::~ChemicalReaction() {
  disconnect();
}

void ChemicalReaction::Initialize() {
  LOG4ESPP_INFO(theLogger, "init ChemicalReaction");
}

void ChemicalReaction::AddReaction(shared_ptr<integrator::Reaction> reaction) {
  reaction->set_dt(dt_);
  reaction->set_interval(interval_);
  reaction->set_rng(rng_);

  // The cutoff of the reaction shouldn't be larger than the cutoff of verletlist.
  if (reaction->cutoff() >  current_cutoff_) reaction->set_cutoff(current_cutoff_);

  reaction_list_.push_back(reaction);
}

void ChemicalReaction::RemoveReaction(int reaction_id) {
  reaction_list_.erase(reaction_list_.begin() + reaction_id);
}

/** Performs all steps of the reactive scheme.
 */
void ChemicalReaction::React() {
  if (integrator->getStep() % (*interval_) != 0) return;

  System& system = getSystemRef();

  LOG4ESPP_INFO(theLogger, "Perform ChemicalReaction");

  *dt_ = integrator->getTimeStep();

  potential_pairs_.clear();
  // loop over VL pairs
  for (PairList::Iterator it(verlet_list_->getPairs()); it.isValid(); ++it) {
    Particle &p1 = *it->first;
    Particle &p2 = *it->second;
    bool found = false;
    int reaction_idx_ = 0;

    for (integrator::ReactionList::iterator it = reaction_list_.begin();
        it != reaction_list_.end() && !found; ++it, ++reaction_idx_) {
      found = (*it)->IsValidPair(p1, p2);
    }

    // If criteria for reaction match, add the indices to potential_pairs_
    if (found) {
      potential_pairs_.insert(std::make_pair(p1.id(), std::make_pair(p2.id(), reaction_idx_)));
    }
  }
  SendMultiMap(potential_pairs_);
  // Here, reduce number of partners to each A to 1
  // Also, keep only non-ghost A
  UniqueA(potential_pairs_);
  SendMultiMap(potential_pairs_);
  // Here, reduce number of partners to each B to 1
  // Also, keep only non-ghost B
  UniqueB(potential_pairs_, effective_pairs_);
  SendMultiMap(effective_pairs_);
  // Use effective_pairs_ to apply the reaction.
  ApplyAR();
}


/** Performs two-way parallel communication to consolidate mm between
neighbours. The parallel scheme is taken from
storage::DomainDecomposition::doGhostCommunication
*/
void ChemicalReaction::SendMultiMap(integrator::ReactionMap &mm) {  //NOLINT
  LOG4ESPP_INFO(theLogger, "Entering sendMultiMap");

  InBuffer inBuffer0(*getSystem()->comm);
  InBuffer inBuffer1(*getSystem()->comm);
  OutBuffer outBuffer(*getSystem()->comm);
  System& system = getSystemRef();
  const storage::NodeGrid& nodeGrid = domdec_->getNodeGrid();

  /* direction loop: x, y, z.
    Here we could in principle build in a one sided ghost
    communication, simply by taking the lr loop only over one
    value. */

  for (int coord = 0; coord < 3; ++coord) {
    /* inverted processing order for ghost force communication,
       since the corner ghosts have to be collected via several
       nodes. We now add back the corner ghost forces first again
       to ghost forces, which only eventually go back to the real
       particle.
    */

    real curCoordBoxL = system.bc->getBoxL()[coord];

    outBuffer.reset();
    // fill outBuffer from mm
    int tmp = mm.size();
    int a, b, c;
    outBuffer.write(tmp);
    for (integrator::ReactionMap::iterator it = mm.begin(); it != mm.end(); it++) {
      a = it->first;
      b = it->second.first;
      c = it->second.second;
      outBuffer.write(a);
      outBuffer.write(b);
      outBuffer.write(c);
    }

    // lr loop: left right
    for (int lr = 0; lr < 2; ++lr) {
      int dir = 2 * coord + lr;
      int oppositeDir = 2 * coord + (1 - lr);
      int dirSize = nodeGrid.getGridSize(coord);
      // Avoids double communication for size 2 directions.
      if ( (dirSize == 2) && (lr == 1) ) continue;

      if (dirSize == 1) {
        LOG4ESPP_DEBUG(theLogger, "no communication");
      } else {
        // prepare send and receive buffers
        longint receiver, sender;
        receiver = nodeGrid.getNodeNeighborIndex(dir);
        sender = nodeGrid.getNodeNeighborIndex(oppositeDir);

        // exchange particles, odd-even rule
        if (nodeGrid.getNodePosition(coord) % 2 == 0) {
          outBuffer.send(receiver, kCrCommTag);
          if (lr == 0)  {
            inBuffer0.recv(sender, kCrCommTag);
          } else {
            inBuffer1.recv(sender, kCrCommTag);
          }
        } else {
          if (lr == 0) {
            inBuffer0.recv(sender, kCrCommTag);
          } else {
            inBuffer1.recv(sender, kCrCommTag);
          }
          outBuffer.send(receiver, kCrCommTag);
        }
      }
    }
    LOG4ESPP_DEBUG(theLogger, "Entering unpack");
    // unpack received data
    // add content of inBuffer to mm
    int lengthA, Aidx, Bidx, Ridx;
    for (int lr = 0; lr < 2; ++lr) {
      int dir = 2 * coord + lr;
      int oppositeDir = 2 * coord + (1 - lr);
      int dirSize = nodeGrid.getGridSize(coord);
      if (dirSize == 1) {
      } else {
      // Avoids double communication for size 2 directions.
        if ( (dirSize == 2) && (lr == 1) ) continue;
        if (lr == 0) {
          inBuffer0.read(lengthA);
        } else {
          inBuffer1.read(lengthA);
        }
        for (longint i = 0; i < lengthA; i++) {
          if (lr == 0) {
            inBuffer0.read(Aidx);
            inBuffer0.read(Bidx);
            inBuffer0.read(Ridx);
          } else {
            inBuffer1.read(Aidx);
            inBuffer1.read(Bidx);
            inBuffer1.read(Ridx);
          }
          mm.insert(std::make_pair(Aidx, std::make_pair(Bidx, Ridx)));
        }
      }
    }
    LOG4ESPP_DEBUG(theLogger, "Leaving unpack");
  }
  LOG4ESPP_INFO(theLogger, "Leaving sendMultiMap");
}

/** Given a multimap mm with several pairs (id1,id2), keep only one pair for
each id1 and return it in place. In addition, only pairs for which
id1 is local are kept.
*/
void ChemicalReaction::UniqueA(integrator::ReactionMap &mm) {  //NOLINT
  // Collect indices
  boost::unordered_set<longint> idxSet;
  boost::unordered_multimap<longint, longint> idxList;
  longint idx1, idx2, reactIdx;
  System& system = getSystemRef();
  boost::unordered_multimap<longint, longint> uniqueList;

  idxSet.clear();
  idxList.clear();
  // Create idxList, containing (idx1, idx2) pairs
  // Create idxSet, containing each idx1 only once
  for (integrator::ReactionMap::iterator it = mm.begin();
          it != mm.end(); it++) {
    idx1 = it->first;
    idx2 = it->second.first;
    reactIdx = it->second.second;
    idxList.insert(std::make_pair(idx1, idx2));
    idxSet.insert(idx1);
  }

  /*
  uniqueList.clear();
  // For each active idx1, pick a partner
  if (idxSet.size() > 0) {
    for (boost::unordered_set<longint>::iterator it = idxSet.begin();
            it != idxSet.end(); it++) {
      idx1 = *it;
      Particle* p = system.storage->lookupLocalParticle(idx1);
      if (p == NULL) continue;
      if (p->ghost()) continue;
      int size = idxList.count(idx1);
      if (size > 0) {
        int pick = (*rng)(size);
        std::pair<
            boost::unordered_multimap<longint, longint>::iterator,
            boost::unordered_multimap<longint, longint>::iterator> candidates;
        candidates = idxList.equal_range(idx1);
        int i = 0;
        for (boost::unordered_multimap<longint, longint>::iterator jt = candidates.first;
                jt != candidates.second; jt++, i++) {
          if (i == pick) {
            uniqueList.insert(std::make_pair(jt->first, jt->second));
            break;
          }
        }
      }
    }
  }
  mm = uniqueList;
  uniqueList.clear();
  */
}

/** Given a multimap mm with several pairs (id1,id2), keep only one pair for
each id2 and return it in place. In addition, only pairs for which
id2 is local are kept.
*/
void ChemicalReaction::UniqueB(
    integrator::ReactionMap& mm, integrator::ReactionMap& nn) {  //NOLINT
  /*
  // Collect indices
  boost::unordered_set<longint> idxSet;
  boost::unordered_multimap<longint, longint> idxList;
  longint idx1, idx2;
  System& system = getSystemRef();
  boost::unordered_multimap<longint, longint> uniqueList;

  idxSet.clear();
  idxList.clear();
  // Create idxList, containing (idx2, idx1) pairs
  // Create idxSet, containing each idx2 only once
  for (boost::unordered_multimap<longint, longint>::iterator it = mm.begin();
          it != mm.end(); it++) {
    idx1 = it->first;
    idx2 = it->second;
    idxList.insert(std::make_pair(idx2, idx1));
    idxSet.insert(idx2);
  }

  uniqueList.clear();
  if (idxSet.size() > 0) {
    // For each active idx1, pick a partner
    for (boost::unordered_set<longint>::iterator it = idxSet.begin();
            it != idxSet.end(); it++) {
      idx2 = *it;
      Particle* p = system.storage->lookupLocalParticle(idx2);
      if (p == NULL) continue;
      if (p->ghost()) continue;
      int size = idxList.count(idx2);
      if (size > 0) {
        int pick = (*rng)(size);
        std::pair<
            boost::unordered_multimap<longint, longint>::iterator,
            boost::unordered_multimap<longint, longint>::iterator> candidates;
        candidates = idxList.equal_range(idx2);
        int i = 0;
        for (boost::unordered_multimap<longint, longint>::iterator jt = candidates.first;
                jt != candidates.second; jt++, i++) {
          if (i == pick) {
            uniqueList.insert(std::make_pair(jt->first, jt->second));
            break;
          }
        }
      }
    }
  }
  nn = uniqueList;
  uniqueList.clear();
  */
}

/** Use the (A,B) list "partners" to add bonds and change the state of the
particles accordingly.
*/
void ChemicalReaction::ApplyAR() {
  /*
  longint A, B, tmp_AB;
  int reaction_idx;
  System& system = getSystemRef();

  LOG4ESPP_INFO(theLogger, "Entering applyAR");

  for (boost::unordered_multimap<longint, longint>::iterator it = effective_pairs_.begin();
          it != effective_pairs_.end(); it++) {
    A = it->second.first;
    reaction_idx = it->second.second;
    B = it->first;

    Reaction reaction = (Reaction)reaction_list_.at(reaction_idx);

    // Change the state of A and B.
    Particle* pA = system.storage->lookupLocalParticle(A);
    Particle* pB = system.storage->lookupLocalParticle(B);
    Particle* tmp;
    if (pA->getType() == reaction.type_b()) {
      tmp = pB;
      tmp_AB = B;
      pB = pA;
      B = A;
      pA = tmp;
      A = tmp_AB;
    }
    if (pA == NULL) {
    } else {
      pA->setState(pA->getState()+reaction.delta_a());
    }
    if (pB == NULL) {
    } else {
      pB->setState(pB->getState()+reaction.delta_b());
    }
    // Add a bond
    fixed_pair_list_->add(A, B);
  }
  LOG4ESPP_INFO(theLogger, "Leaving applyAR");
  */
}


void ChemicalReaction::disconnect() {
  initialize_.disconnect();
  react_.disconnect();
}

void ChemicalReaction::connect() {
  // connect to initialization inside run()
  initialize_ = integrator->runInit.connect(
    boost::bind(&ChemicalReaction::Initialize, this));

  react_ = integrator->aftIntV.connect(
    boost::bind(&ChemicalReaction::React, this));
}



/****************************************************
 ** REGISTRATION WITH PYTHON
 ****************************************************/
void ChemicalReaction::registerPython() {
    using namespace espresso::python;  //NOLINT
    class_<ChemicalReaction, shared_ptr<ChemicalReaction>, bases<Extension> >
      ("integrator_ChemicalReaction",
          init<shared_ptr<System>,
          shared_ptr<VerletList>,
          shared_ptr<FixedPairList>,
          shared_ptr<storage::DomainDecomposition> >())
        .def("connect", &ChemicalReaction::connect)
        .def("disconnect", &ChemicalReaction::disconnect)
        .def("addReaction", &ChemicalReaction::AddReaction)
        .def("removeReaction", &ChemicalReaction::RemoveReaction)
        .add_property(
          "interval",
          &ChemicalReaction::interval,
          &ChemicalReaction::set_interval);
    }
}  // namespace integrator
}  // namespace espresso

