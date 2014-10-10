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

#include <cstdio>
#include <utility>
#include <vector>
#include <set>

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

/** Checks if the particles pair is valid. */
bool Reaction::IsValidPair(const Particle& p1, const Particle& p2) {
  if (IsValidState(p1, p2)) {
    Real3D distance = p1.position() - p2.position();
    real distance_2 = distance.sqr();

    return ((distance_2 < cutoff_sqr_) && ((*rng_)() < rate_*(*dt_)*(*interval_)));
  }
  return false;
}


/** Checks if the particles has correct state. */
bool Reaction::IsValidState(const Particle& p1, const Particle& p2) {
  if ((p1.res_id() == p2.res_id()) && !intramolecular_)
    return false;

  int p1_state = p1.state();
  int p2_state = p2.state();
  if ((p1.type() == type_a_) && (p2.type() == type_b_) &&
      (p1_state >= min_state_a_ && p1_state < max_state_a_) &&
      (p2_state >= min_state_b_ && p2_state < max_state_b_)) {
    return true;
  } else if ((p1.type() == type_b_) && (p2.type() == type_a_) &&
             (p1_state >= min_state_b_ && p1_state < max_state_b_) &&
             (p2_state >= min_state_a_ && p2_state < max_state_a_)) {
    return true;
  }
  return false;
}


void Reaction::registerPython() {
  using namespace espresso::python;  //NOLINT
  class_< Reaction, shared_ptr<Reaction> >
    ("integrator_Reaction",
     init<int, int, int, int, int, int, int, int, real, real, bool>())
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
        &Reaction::max_state_a,
        &Reaction::set_max_state_a)
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
        &Reaction::set_cutoff)
      .add_property(
        "intramolecular",
        &Reaction::intramolecular,
        &Reaction::set_intramolecular);
}

/*
void SynthesisReaction::registerPython() {
  using namespace espresso::python;  //NOLINT
  class_< SynthesisReaction, bases<Reaction>, shared_ptr<Reaction> >
    ("integrator_SynthesisReaction",
     init<int, int, int, int, int, int, int, int, real, real, bool>());
}
*/


/**
 * ChemicalReaction
 * */
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
  effective_pairs_.clear();
  // loop over VL pairs
  for (PairList::Iterator it(verlet_list_->getPairs()); it.isValid(); ++it) {
    Particle &p1 = *it->first;
    Particle &p2 = *it->second;
    int reaction_idx_ = 0;

    for (integrator::ReactionList::iterator it = reaction_list_.begin();
        it != reaction_list_.end(); ++it, ++reaction_idx_) {
      if ((*it)->IsValidPair(p1, p2)) {
        potential_pairs_.insert(std::make_pair(p1.id(), std::make_pair(p2.id(), reaction_idx_)));
      }
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

  // Prepare out buffer with the reactions that potential will happen on this node.
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

  /* direction loop: x, y, z.
    Here we could in principle build in a one sided ghost
    communication, simply by taking the lr loop only over one
    value. */

  for (int direction = 0; direction < 3; ++direction) {
    /* inverted processing order for ghost force communication,
       since the corner ghosts have to be collected via several
       nodes. We now add back the corner ghost forces first again
       to ghost forces, which only eventually go back to the real
       particle.
    */
    real curCoordBoxL = system.bc->getBoxL()[direction];

    // lr loop: left right
    for (int left_right_dir = 0; left_right_dir < 2; ++left_right_dir) {
      int neighbour_node = 2 * direction + left_right_dir;
      int opposite_neighbour_node = 2 * direction + (1 - left_right_dir);
      int direction_size = nodeGrid.getGridSize(direction);

      // Avoids double communication for size 2 directions.
      if ( (direction_size == 2) && (left_right_dir == 1) ) continue;

      if (direction_size == 1) {  // No communication needed in this direction.
        LOG4ESPP_DEBUG(theLogger, "No communication needed.");
      } else {
        // prepare send and receive buffers
        longint receiver, sender;
        receiver = nodeGrid.getNodeNeighborIndex(neighbour_node);
        sender = nodeGrid.getNodeNeighborIndex(opposite_neighbour_node);

        // exchange particles, odd-even rule
        if (nodeGrid.getNodePosition(direction) % 2 == 0) {
          outBuffer.send(receiver, kCrCommTag);
          if (left_right_dir == 0)  {
            inBuffer0.recv(sender, kCrCommTag);
          } else {
            inBuffer1.recv(sender, kCrCommTag);
          }
        } else {
          if (left_right_dir == 0) {
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
    int data_length, idx_a, idx_b, reaction_idx;

    for (int left_right_dir = 0; left_right_dir < 2; ++left_right_dir) {
      int neighbour_node = 2 * direction + left_right_dir;
      int opposite_neighbour_node = 2 * direction + (1 - left_right_dir);
      int direction_size = nodeGrid.getGridSize(direction);
      if (direction_size > 1) {
        // Avoids double communication for size 2 directions.
        if ( (direction_size == 2) && (left_right_dir == 1) ) continue;
        if (left_right_dir == 0) {
          inBuffer0.read(data_length);
        } else {
          inBuffer1.read(data_length);
        }
        for (longint i = 0; i < data_length; i++) {
          if (left_right_dir == 0) {
            inBuffer0.read(idx_a);
            inBuffer0.read(idx_b);
            inBuffer0.read(reaction_idx);
          } else {
            inBuffer1.read(idx_a);
            inBuffer1.read(idx_b);
            inBuffer1.read(reaction_idx);
          }
          mm.insert(std::make_pair(idx_a, std::make_pair(idx_b, reaction_idx)));
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
void ChemicalReaction::UniqueA(integrator::ReactionMap &potential_candidates) { //NOLINT
  LOG4ESPP_INFO(theLogger, "UniqueA");
  System& system = getSystemRef();
  integrator::ReactionMap unique_list_of_candidates;
  boost::unordered_set<longint> a_indexes;

  unique_list_of_candidates.clear();
  // Gets the list of indexes of particle a. REQOPT
  for (integrator::ReactionMap::iterator it = potential_candidates.begin();
       it != potential_candidates.end(); ++it) {
    a_indexes.insert(it->first);
  }

  // For each active idx1, pick a partner
  if (a_indexes.size() > 0) {
    int idx_a;
    real max_reaction_rate;
    // rate => idx_b, reaction_id
    boost::unordered_multimap<real, std::pair<longint, int> > rate_idx_b;
    boost::unordered_multimap<real, std::pair<longint, int> >::local_iterator idx_b_reaction_id;
    Particle *p = NULL;
    std::pair<integrator::ReactionMap::iterator,
              integrator::ReactionMap::iterator> candidates_b;
    // Iterate over the ids of particle A, looking for the particle B
    for (boost::unordered_set<longint>::iterator it = a_indexes.begin();
              it != a_indexes.end(); ++it) {
      idx_a = *it;
      p = system.storage->lookupLocalParticle(idx_a);

      if (p == NULL || p->ghost()) continue;

      // Select all possible candidates (so
      candidates_b = potential_candidates.equal_range(idx_a);

      // Group the candidates by the reaction rate.
      max_reaction_rate = -1;
      rate_idx_b.clear();
      for (integrator::ReactionMap::iterator jt = candidates_b.first;
              jt != candidates_b.second; ++jt) {
        integrator::Reaction reaction = *reaction_list_.at(jt->second.second);
        real reaction_rate = reaction.rate();
        if (reaction_rate > max_reaction_rate) {
          max_reaction_rate = reaction_rate;
        }
        // rate => (idx_b, reaction_id)
        rate_idx_b.insert(
            std::make_pair(reaction_rate, std::make_pair(jt->second.first, jt->second.second)));
      }

      // Found reaction with the maximum rate. If there are several candidates with the same
      // rate, then we choose randomly.
      if (max_reaction_rate != -1) {
        int bucket_size = rate_idx_b.count(max_reaction_rate);
        // Pick up random number in given range.
        int pick_offset = (*rng_)(bucket_size);
        idx_b_reaction_id = rate_idx_b.begin(rate_idx_b.bucket(max_reaction_rate));
        std::advance(idx_b_reaction_id, pick_offset);

        unique_list_of_candidates.insert(std::make_pair(idx_a, idx_b_reaction_id->second));
      }
    }
  }
  potential_candidates = unique_list_of_candidates;
}

/** Given a multimap mm with several pairs (id1,id2), keep only one pair for
each id2 and return it in place. In addition, only pairs for which
id2 is local are kept.
*/
void ChemicalReaction::UniqueB(
    integrator::ReactionMap& potential_candidates,  //NOLINT
    integrator::ReactionMap& effective_candidates) {  //NOLINT

  LOG4ESPP_INFO(theLogger, "UniqueB");

  typedef boost::unordered_set<longint> Indexes;
  typedef boost::unordered_multimap<real, std::pair<longint, int> > RateParticleIdx;

  System& system = getSystemRef();
  Indexes b_indexes;
  integrator::ReactionMap reverse_candidates;

  effective_candidates.clear();

  // Collect the b particle pairs. REQOPT
  for (integrator::ReactionMap::iterator it = potential_candidates.begin();
       it != potential_candidates.end(); ++it) {
    b_indexes.insert(it->second.first);
    reverse_candidates.insert(
        std::make_pair(
          it->second.first,
          std::make_pair(it->first, it->second.second)));
  }

  if (b_indexes.size() > 0) {
    int idx_b;
    real max_reaction_rate;
    // rate => idx_a, reaction_id
    RateParticleIdx rate_idx_a;
    RateParticleIdx::local_iterator idx_a_reaction_id;
    std::pair<integrator::ReactionMap::iterator,
              integrator::ReactionMap::iterator> candidates_a;
    Particle *p = NULL;

    for (Indexes::iterator it = b_indexes.begin(); it != b_indexes.end(); ++it) {
      idx_b = *it;
      p = system.storage->lookupLocalParticle(idx_b);

      if (p == NULL || p->ghost()) continue;

      candidates_a = reverse_candidates.equal_range(idx_b);
      max_reaction_rate = -1;
      rate_idx_a.clear();
      for (integrator::ReactionMap::iterator jt = candidates_a.first;
          jt != candidates_a.second; ++jt) {
        integrator::Reaction reaction = *reaction_list_.at(jt->second.second);
        real reaction_rate = reaction.rate();
        if (reaction_rate > max_reaction_rate) {
          max_reaction_rate = reaction_rate;
        }
        rate_idx_a.insert(
            std::make_pair(reaction_rate, std::make_pair(jt->second.first, jt->second.second)));
      }
      // Found reaction with the maximum rate. If there are several candidates
      // then select randomly.
      if (max_reaction_rate > -1) {
        int bucket_size = rate_idx_a.count(max_reaction_rate);
        int pick_offset = (*rng_)(bucket_size);
        idx_a_reaction_id = rate_idx_a.begin(rate_idx_a.bucket(max_reaction_rate));
        std::advance(idx_a_reaction_id, pick_offset);
        effective_candidates.insert(
            std::make_pair(
              idx_a_reaction_id->second.first,
              std::make_pair(idx_b, idx_a_reaction_id->second.second)));
      }
    }
  }
}


/** Use the (A,B) list "partners" to add bonds and change the state of the
particles accordingly.
*/
void ChemicalReaction::ApplyAR() {
  Particle* pA = NULL;
  Particle* pB = NULL;
  System& system = getSystemRef();

  LOG4ESPP_INFO(theLogger, "Entering applyAR");

  for (integrator::ReactionMap::iterator it = effective_pairs_.begin();
          it != effective_pairs_.end(); it++) {
    integrator::Reaction reaction = *reaction_list_.at(it->second.second);

    // Change the state of A and B.
    pA = system.storage->lookupLocalParticle(it->first);
    pB = system.storage->lookupLocalParticle(it->second.first);
    if (pA != NULL && pB != NULL) {
      if (reaction.IsValidState(*pA, *pB)) {
        if (pA->getType() == reaction.type_a()) {
          pA->setState(pA->getState() + reaction.delta_a());
          pB->setState(pB->getState() + reaction.delta_b());
          pB->setResId(pA->getResId());
        } else if (pA->getType() == reaction.type_b()) {
          pA->setState(pA->getState() + reaction.delta_b());
          pB->setState(pB->getState() + reaction.delta_a());
          pA->setResId(pB->getResId());
        }
        fixed_pair_list_->add(it->first, it->second.first);  // The order does not matter.
        LOG4ESPP_DEBUG(theLogger, "Created pair.");
      }
    }
  }
  LOG4ESPP_INFO(theLogger, "Leaving applyAR");
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

