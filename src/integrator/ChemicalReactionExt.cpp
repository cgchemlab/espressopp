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

#include "ChemicalReactionExt.hpp"

#include <utility>
#include <set>
#include <numeric>
#include <math.h>

#include "storage/Storage.hpp"
#include "iterator/CellListIterator.hpp"
#include "esutil/RNG.hpp"
#include "bc/BC.hpp"
#include "storage/NodeGrid.hpp"
#include "storage/DomainDecomposition.hpp"
#include "FixDistances.hpp"

namespace espressopp {
namespace integrator {
LOG4ESPP_LOGGER(ChemicalReaction::theLogger, "ChemicalReaction");

/** ChemicalReaction extension to the integrator
 *
 * @param system The espressopp.System object.
 * @param verletList The espressopp.VerletList object.
 * @param domdec The espressopp.storage.DomainDecomposition object.
 */
ChemicalReaction::ChemicalReaction(shared_ptr<System> system, shared_ptr<VerletList> verletList,
    shared_ptr<storage::DomainDecomposition> domdec, shared_ptr<TopologyManager> tm)
      :Extension(system),
          verlet_list_(verletList),
          domdec_(domdec), tm_(tm), is_nearest_(false) {
  type = Extension::Reaction;

  current_cutoff_ = verletList->getVerletCutoff() - system->getSkin();

  if (!system->rng)
    throw std::runtime_error("System has no RNG.");

  rng_ = system->rng;
  LOG4ESPP_INFO(theLogger, "ChemicalReaction constructed");
  dt_ = boost::make_shared<real>();
  interval_ = boost::make_shared<int>();

  reaction_list_ = ReactionList();
  reverse_reaction_list_ = ReactionList();

  bond_limit_ = -1;

  resetTimers();
}

ChemicalReaction::~ChemicalReaction() {
  LOG4ESPP_INFO(theLogger, "Destructor ChemicalReaction");
  disconnect();
}

/** Adds the chemical reaction to the list of reactions
 *
 * @param reaction espressopp.Reaction object.
 */
void ChemicalReaction::addReaction(boost::shared_ptr<integrator::Reaction> reaction) {
  if (!reaction->reaction_cutoff())
    throw std::runtime_error("Reaction object does not have ReactionCutoff object.");

  reaction->set_dt(dt_);
  reaction->set_interval(interval_);
  reaction->set_rng(rng_);
  reaction->setTopologyManager(tm_);

  bc::BC &bc = *getSystemRef().bc;

  reaction->set_bc(&bc);

  if (!reaction->reverse()) {
    // If VL cutoff is smaller than reaction, increase it.
    if (reaction->cutoff() > current_cutoff_) {
      LOG4ESPP_INFO(theLogger, "VL cutoff is extended to match with reaction cutoff");
      verlet_list_->setVerletCutoff(reaction->cutoff());
      current_cutoff_ = reaction->cutoff();
    }

    LOG4ESPP_INFO(theLogger, "Added reaction");
    reaction_list_.push_back(reaction);
  } else {
    LOG4ESPP_INFO(theLogger, "Add reverse reaction");

    // In this case, VL cutoff does not matter. Adds reaction on separate list.
    reverse_reaction_list_.push_back(reaction);
  }
}

/** Performs all steps of the reactive scheme. */
void ChemicalReaction::React() {
  if (integrator->getStep() % (*interval_) != 0)
    return;

  System &system = getSystemRef();

  LOG4ESPP_DEBUG(theLogger, "Perform ChemicalReaction");

  *dt_ = integrator->getTimeStep();

  potential_pairs_.clear();
  effective_pairs_.clear();

  wallTimer.startMeasure();
  // loop over VL pairs
  for (PairList::Iterator it(verlet_list_->getPairs()); it.isValid(); ++it) {
    Particle &p1 = *it->first;
    Particle &p2 = *it->second;
    int reaction_idx_ = 0;

    for (ReactionList::iterator it = reaction_list_.begin(); it != reaction_list_.end(); ++it, ++reaction_idx_) {
      if (!(*it)->active())  // if raction is not active, skip it.
        continue;

      ReactedPair p;

      if ((*it)->isValidPair(p1, p2, p)) {
        potential_pairs_.insert(
            std::make_pair(p.first->id(),
                           std::make_pair(
                               p.second->id(),
                               ReactionDef(reaction_idx_, p.reaction_rate, p.r_sqr)
                           )
            )
        );
      }
    }
  } // end loop over VL pairs
  timeLoopPair += wallTimer.stopMeasure();

  wallTimer.startMeasure();
  sendMultiMap(potential_pairs_);

  // Here, reduce number of partners to each A to 1
  // Also, keep only non-ghost A
  UniqueA(potential_pairs_);
  sendMultiMap(potential_pairs_);
  // Here, reduce number of partners to each B to 1
  // Also, keep only non-ghost B
  UniqueB(potential_pairs_, effective_pairs_);
  // Distribute effective pairs
  sendMultiMap(effective_pairs_);

  // Use effective_pairs_ to apply the reaction.
  std::set<Particle *> modified_particles;
  timeComm += wallTimer.stopMeasure();

  wallTimer.startMeasure();
  // First, remove pairs.
  ApplyDR(modified_particles);
  timeApplyDR += wallTimer.stopMeasure();

  // Synchronize, all cpus should finish dissocition part.
  (*system.comm).barrier();

  wallTimer.startMeasure();
  // Now, accept new pairs.
  ApplyAR(modified_particles);
  timeApplyAR += wallTimer.stopMeasure();

  // Synchronize, all cpus should finish association part.
  (*system.comm).barrier();

  wallTimer.startMeasure();
  // Update the ghost particles.
  updateGhost(modified_particles);
  timeUpdateGhost += wallTimer.stopMeasure();

  LOG4ESPP_DEBUG(theLogger, "Finished react()");
  LOG4ESPP_DEBUG(theLogger, "Leaving react()");
}

/** Performs two-way parallel communication to consolidate mm between
   neighbours. The parallel scheme is taken from
   storage::DomainDecomposition::doGhostCommunication
 */

void ChemicalReaction::printMultiMap(ReactionMap &rmap, std::string comment) {
  System &system = getSystemRef();
  for (integrator::ReactionMap::iterator it = rmap.begin(); it != rmap.end(); it++) {
    std::cout << comment << "mm on\t" << system.comm->rank() << "\t" << it->first << "\t" << it->second.first
        << "\t" << it->second.second.reaction_id << "\t"
        << it->second.second.reaction_rate << "\t" << it->second.second.reaction_r_sqr << std::endl;
  }
}

void ChemicalReaction::sendMultiMap(integrator::ReactionMap &mm) {// NOLINT
  LOG4ESPP_DEBUG(theLogger, "Entering sendMultiMap");

  System &system = getSystemRef();

  InBuffer in_buffer_0(*system.comm);
  InBuffer in_buffer_1(*system.comm);
  OutBuffer out_buffer(*system.comm);
  const storage::NodeGrid &node_grid = domdec_->getNodeGrid();

  // Prepare out buffer with the reactions that potential will happen on this node.
  out_buffer.reset();
  in_buffer_0.reset();
  in_buffer_1.reset();

  // Fill out_buffer from mm.
  int array_size = mm.size();
  int a, b, c;

  real d, r_sqr;

  out_buffer.write(array_size);

  for (integrator::ReactionMap::iterator it = mm.begin(); it != mm.end();
      it++) {
    a = it->first;  // particle id
    b = it->second.first;  // particle id
    c = it->second.second.reaction_id;  // reaction id
    d = it->second.second.reaction_rate; // reaction rate for this pair.
    r_sqr = it->second.second.reaction_r_sqr;  // reaction distance for this pair.
    out_buffer.write(a);
    out_buffer.write(b);
    out_buffer.write(c);
    out_buffer.write(d);
    out_buffer.write(r_sqr);

  }

  LOG4ESPP_DEBUG(theLogger, "OutBuffer.size=" << out_buffer.getSize());

  /* direction loop: x, y, z.
     Here we could in principle build in a one sided ghost
     communication, simply by taking the lr loop only over one
     value. */

  int data_length, idx_a, idx_b, reaction_idx, direction_size;
  real reaction_rate, reaction_r_sqr;

  for (int direction = 0; direction < 3; ++direction) {
    /* inverted processing order for ghost force communication,
       since the corner ghosts have to be collected via several
       nodes. We now add back the corner ghost forces first again
       to ghost forces, which only eventually go back to the real
       particle.
     */
    direction_size = node_grid.getGridSize(direction);

    if (direction_size == 1) {
      LOG4ESPP_DEBUG(theLogger, "No communication needed.");
      continue;
    }

    // lr loop: left right
    for (int left_right_dir = 0; left_right_dir < 2; ++left_right_dir) {
      // Avoids double communication for size 2 directions.
      if ((direction_size == 2) && (left_right_dir == 1))
        continue;

      // prepare send and receive buffers
      longint receiver, sender;

      receiver = node_grid.getNodeNeighborIndex(2 * direction + left_right_dir);
      sender = node_grid.getNodeNeighborIndex(2 * direction + (1 - left_right_dir));

      // exchange particles, odd-even rule. getNodePosition returns the position
      // of the current node.
      if (node_grid.getNodePosition(direction) % 2 == 0) {
        // sending the data
        out_buffer.send(receiver, kCrCommTag);

        // receiving the data
        if (left_right_dir == 0) {
          in_buffer_0.recv(sender, kCrCommTag);
        } else {
          in_buffer_1.recv(sender, kCrCommTag);
        }
      } else {
        // receiving the data
        if (left_right_dir == 0) {
          in_buffer_0.recv(sender, kCrCommTag);
        } else {
          in_buffer_1.recv(sender, kCrCommTag);
        }

        // sending the data
        out_buffer.send(receiver, kCrCommTag);
      }
    }

    LOG4ESPP_DEBUG(theLogger, "Entering unpack");

    // Unpacking phase. Get the content from buffers and put in output multimap.
    for (int left_right_dir = 0; left_right_dir < 2; ++left_right_dir) {
      // Avoids double communication for size 2 directions.
      if ((direction_size == 2) && (left_right_dir == 1))
        continue;

      if (left_right_dir == 0) {
        in_buffer_0.read(data_length);
      } else {
        in_buffer_1.read(data_length);
      }

      for (longint i = 0; i < data_length; i++) {
        if (left_right_dir == 0) {
          in_buffer_0.read(idx_a);
          in_buffer_0.read(idx_b);
          in_buffer_0.read(reaction_idx);
          in_buffer_0.read(reaction_rate);
          in_buffer_0.read(reaction_r_sqr);
        } else {
          in_buffer_1.read(idx_a);
          in_buffer_1.read(idx_b);
          in_buffer_1.read(reaction_idx);
          in_buffer_1.read(reaction_rate);
          in_buffer_1.read(reaction_r_sqr);
        }

        mm.insert(
            std::make_pair(
                idx_a,
                std::make_pair(
                    idx_b,
                    ReactionDef(reaction_idx, reaction_rate, reaction_r_sqr)
                )
            )
        );
      }
    }

    LOG4ESPP_DEBUG(theLogger, "Leaving unpack");
  }

  LOG4ESPP_DEBUG(theLogger, "Leaving sendMultiMap");
}

/** Performs two-way parallel communication to update the ghost particles.
 * The parallel scheme is taken from
 * storage::DomainDecomposition::doGhostCommunication
 */
void ChemicalReaction::updateGhost(const std::set<Particle *> &modified_particles) {// NOLINT
  LOG4ESPP_DEBUG(theLogger, "Entering updateGhost");

  System &system = getSystemRef();

  InBuffer in_buffer_0(*getSystem()->comm);
  InBuffer in_buffer_1(*getSystem()->comm);
  OutBuffer out_buffer(*getSystem()->comm);
  const storage::NodeGrid &node_grid = domdec_->getNodeGrid();

  // Prepare out buffer with the reactions that potential will happen on this node.
  out_buffer.reset();
  in_buffer_0.reset();
  in_buffer_1.reset();

  // Fill out_buffer from the particles properties.
  longint data_length = modified_particles.size();
  longint p_id, p_type, p_res_id, p_state;
  real p_mass, p_q, p_lambda;

  out_buffer.write(data_length);

  for (std::set<Particle *>::const_iterator it = modified_particles.begin();
      it != modified_particles.end();
      ++it) {
    p_id = (*it)->id();
    p_type = (*it)->type();
    p_mass = (*it)->mass();
    p_q = (*it)->q();
    p_res_id = (*it)->res_id();
    p_lambda = (*it)->lambda();
    p_state = (*it)->state();
    out_buffer.write(p_id);
    out_buffer.write(p_type);
    out_buffer.write(p_mass);
    out_buffer.write(p_q);
    out_buffer.write(p_res_id);
    out_buffer.write(p_lambda);
    out_buffer.write(p_state);
  }

  // Temporary data.
  Particle *particle = NULL;
  int direction_size = 0;

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
    direction_size = node_grid.getGridSize(direction);

    if (direction_size == 1) {
      LOG4ESPP_DEBUG(theLogger, "No communication needed.");
      continue;
    }

    // lr loop: left right
    for (int left_right_dir = 0; left_right_dir < 2; ++left_right_dir) {
      // Avoids double communication for size 2 directions.
      if ((direction_size == 2) && (left_right_dir == 1))
        continue;

      // prepare send and receive buffers
      longint receiver, sender;

      receiver = node_grid.getNodeNeighborIndex(2 * direction + left_right_dir);
      sender = node_grid.getNodeNeighborIndex(2 * direction + (1 - left_right_dir));

      // exchange particles, odd-even rule. getNodePosition returns the position
      // of the current node.
      if (node_grid.getNodePosition(direction) % 2 == 0) {
        // sending the data
        out_buffer.send(receiver, kCrCommTag);

        // receiving the data
        if (left_right_dir == 0) {
          in_buffer_0.recv(sender, kCrCommTag);
        } else {
          in_buffer_1.recv(sender, kCrCommTag);
        }
      } else {
        // receiving the data
        if (left_right_dir == 0) {
          in_buffer_0.recv(sender, kCrCommTag);
        } else {
          in_buffer_1.recv(sender, kCrCommTag);
        }

        // sending the data
        out_buffer.send(receiver, kCrCommTag);
      }
    }

    LOG4ESPP_DEBUG(theLogger, "Entering unpack");

    // Unpacking phase. Get the parameters and set the data to particles.
    for (int left_right_dir = 0; left_right_dir < 2; ++left_right_dir) {
      // Avoids double communication for size 2 directions.
      if ((direction_size == 2) && (left_right_dir == 1))
        continue;

      if (left_right_dir == 0) {
        in_buffer_0.read(data_length);
      } else {
        in_buffer_1.read(data_length);
      }

      for (longint i = 0; i < data_length; i++) {
        if (left_right_dir == 0) {
          in_buffer_0.read(p_id);
          in_buffer_0.read(p_type);
          in_buffer_0.read(p_mass);
          in_buffer_0.read(p_q);
          in_buffer_0.read(p_res_id);
          in_buffer_0.read(p_lambda);
          in_buffer_0.read(p_state);
        } else {
          in_buffer_1.read(p_id);
          in_buffer_1.read(p_type);
          in_buffer_1.read(p_mass);
          in_buffer_1.read(p_q);
          in_buffer_1.read(p_res_id);
          in_buffer_1.read(p_lambda);
          in_buffer_1.read(p_state);
        }

        // Update the ghost particle data on neighbour CPUs.
        particle = system.storage->lookupLocalParticle(p_id);

        if (particle != NULL) {
          LOG4ESPP_DEBUG(theLogger, "Update particle data");
          particle->setType(p_type);
          particle->setMass(p_mass);
          particle->setQ(p_q);
          particle->setResId(p_res_id);
          particle->setLambda(p_lambda);
          particle->setState(p_state);
        }
      }
    }

    LOG4ESPP_DEBUG(theLogger, "Leaving unpack");
  }

  LOG4ESPP_DEBUG(theLogger, "Leaving updateGhost");
}

/** Given a multimap mm with several pairs (id1,id2), keep only one pair for
   each id1 and return it in place. In addition, only pairs for which
   id1 is local are kept.
 */
void ChemicalReaction::UniqueA(integrator::ReactionMap &potential_candidates) {// NOLINT
  LOG4ESPP_DEBUG(theLogger, "UniqueA");

  System &system = getSystemRef();
  integrator::ReactionMap unique_list_of_candidates;
  boost::unordered_set<longint> a_indexes;

  Particle *p = NULL;

  unique_list_of_candidates.clear();

  // Gets the list of indexes of particle a. Gets only real particles, skip ghost.
  for (integrator::ReactionMap::iterator it = potential_candidates.begin();
      it != potential_candidates.end(); ++it) {
    p = system.storage->lookupRealParticle(it->first);

    if (p == NULL)
      continue;

    a_indexes.insert(it->first);
  }

  // For each active idx1, pick a partner
  if (a_indexes.size() > 0) {
    int idx_a;
    real rc;
    real max_rc;

    // reaction_coordinate => idx_b, ReactionDef(reaction_id, reaction_rate) or r_sqr distance
    typedef boost::unordered_multimap<real, std::pair<longint, ReactionDef > > LocalRateIdx;
    LocalRateIdx rc_idx_b;
    LocalRateIdx::local_iterator idx_b_reaction_id;

    // Iterators for the equal_range.
    std::pair<integrator::ReactionMap::iterator,
              integrator::ReactionMap::iterator> candidates_b;

    // Iterate over the ids of particle A, looking for the particle B
    for (boost::unordered_set<longint>::iterator it = a_indexes.begin();
        it != a_indexes.end(); ++it) {
      idx_a = *it;

      // Group the candidates by the reaction rate.
      if (is_nearest_)  // silly
        max_rc = 10e18;
      else
        max_rc = -1;

      // Select all possible candidates
      candidates_b = potential_candidates.equal_range(idx_a);

      rc_idx_b.clear();

      for (integrator::ReactionMap::iterator jt = candidates_b.first;
          jt != candidates_b.second; ++jt) {
        boost::shared_ptr<integrator::Reaction> reaction = reaction_list_.at(jt->second.second.reaction_id);
        real reaction_rate = jt->second.second.reaction_rate;
        real reaction_r_sqr = jt->second.second.reaction_r_sqr;

        if (is_nearest_) {
          if (reaction_r_sqr < max_rc)
            max_rc = reaction_r_sqr;
        } else {
          if (reaction_rate > max_rc)
            max_rc = reaction_rate;
        }

        // Use reaction coordinate, distance of reaction rate.
        if (is_nearest_)
          rc = reaction_r_sqr;
        else
          rc = reaction_rate;

        // rc => (idx_b, reaction_id)
        rc_idx_b.insert(std::make_pair(rc, std::make_pair(jt->second.first, jt->second.second)));
      }

      // Found reaction with the maximum rate. If there are several candidates with the same
      // rate, then we choose randomly.
      if (max_rc != -1) {
        int bucket_size = rc_idx_b.count(max_rc);

        int pick_offset = 0;
        if (bucket_size > 1)
          // Pick up random number in given range.
          pick_offset = (*rng_)(bucket_size);

        idx_b_reaction_id = rc_idx_b.begin(rc_idx_b.bucket(max_rc));

        std::advance(idx_b_reaction_id, pick_offset);

        unique_list_of_candidates.insert(
            std::make_pair(idx_a,
                           std::make_pair(
                               idx_b_reaction_id->second.first,
                               idx_b_reaction_id->second.second
                           )));
      }
    }
  }

  //@todo(jakub): I'm not sure if this is an efficient approach.
  potential_candidates = unique_list_of_candidates;
}

/** Given a multimap mm with several pairs (id1,id2), keep only one pair for
   each id2 and return it in place. In addition, only pairs for which
   id2 is local are kept.
 */
void ChemicalReaction::UniqueB(integrator::ReactionMap &potential_candidates,// NOLINT
                               integrator::ReactionMap &effective_candidates) {// NOLINT
  LOG4ESPP_DEBUG(theLogger, "UniqueB");

  typedef boost::unordered_set<longint> Indexes;
  typedef boost::unordered_multimap<real, std::pair<longint, ReactionDef> > RateParticleIdx;

  System &system = getSystemRef();
  Indexes b_indexes;
  integrator::ReactionMap reverse_candidates;

  Particle *p = NULL;

  effective_candidates.clear();

  // Collect the b particle pairs. REQOPT
  for (integrator::ReactionMap::iterator it = potential_candidates.begin();
      it != potential_candidates.end(); ++it) {
    p = system.storage->lookupRealParticle(it->second.first);

    if (p == NULL)
      continue;

    b_indexes.insert(it->second.first);
    reverse_candidates.insert(
        std::make_pair(it->second.first,  // idB
                       std::make_pair(
                           it->first,     // idA
                           it->second.second)));  // ReactionDef
  }

  if (b_indexes.size() > 0) {
    int idx_b;
    real max_rc;
    real rc;

    // rate => idx_a, reaction_id, reaction_rate
    RateParticleIdx rc_idx_a;
    RateParticleIdx::local_iterator idx_a_reaction_id;
    std::pair<integrator::ReactionMap::iterator,
    integrator::ReactionMap::iterator> candidates_a;

    for (Indexes::iterator it = b_indexes.begin(); it != b_indexes.end(); ++it) {
      idx_b = *it;

      candidates_a = reverse_candidates.equal_range(idx_b);
      // Group the candidates by the reaction rate.
      if (is_nearest_)  // silly
        max_rc = 10e18;
      else
        max_rc = -1;

      rc_idx_a.clear();

      for (integrator::ReactionMap::iterator jt = candidates_a.first;
          jt != candidates_a.second; ++jt) {
        boost::shared_ptr<integrator::Reaction> reaction = reaction_list_.at(jt->second.second.reaction_id);
        real reaction_rate = jt->second.second.reaction_rate;
        real reaction_r_sqr = jt->second.second.reaction_r_sqr;

        if (is_nearest_) {
          if (reaction_r_sqr < max_rc)
            max_rc = reaction_r_sqr;
        } else {
          if (reaction_rate > max_rc)
            max_rc = reaction_rate;
        }

        // Use reaction coordinate, distance of reaction rate.
        if (is_nearest_)
          rc = reaction_r_sqr;
        else
          rc = reaction_rate;

        rc_idx_a.insert(
            std::make_pair(rc, std::make_pair(jt->second.first, jt->second.second)));
      }

      // Found reaction with the maximum rate. If there are several candidates
      // then select randomly.
      if (max_rc != -1) {
        int bucket_size = rc_idx_a.count(max_rc);
        int pick_offset = 0;
        if (bucket_size > 1)
          pick_offset = (*rng_)(bucket_size);

        idx_a_reaction_id = rc_idx_a.begin(rc_idx_a.bucket(max_rc));

        std::advance(idx_a_reaction_id, pick_offset);

        effective_candidates.insert(
            std::make_pair(
                idx_a_reaction_id->second.first,
                std::make_pair(idx_b,
                               idx_a_reaction_id->second.second)));
      }
    }
  }
}

/** ApplyDR to remove bonds and change the state of the particles
 * accordingly
 */
void ChemicalReaction::ApplyDR(std::set<Particle *> &modified_particles) {
  LOG4ESPP_DEBUG(theLogger, "Entering applyDR");

  // Iterate over reverse reaction. For every reaction, iterate overy particles pairs and
  // decide to remove or keep bond.
  std::set<Particle *> tmp;

  for (ReactionList::iterator it = reverse_reaction_list_.begin();
      it != reverse_reaction_list_.end(); ++it) {
    integrator::Reaction &r = **it;

    if (!r.active())
      continue;

    bool updated_fpl = false;

    for (FixedPairList::PairList::Iterator itp(*((*it)->fixed_pair_list_)); itp.isValid(); ++itp) {
      Particle &p1 = *itp->first;
      Particle &p2 = *itp->second;
      ReactedPair p;

      if ((*it)->isValidPair(p1, p2, p)) {
        // Ok we will remove this pair.
        if (!p1.ghost() && !p1.ghost())
          (*it)->fixed_pair_list_->remove(p1.id(), p2.id());
        else if (p2.ghost())
          (*it)->fixed_pair_list_->remove(p1.id(), p2.id());
        else
          (*it)->fixed_pair_list_->remove(p2.id(), p1.id());

        // Run postprocess.
        p.first->setState(p.first->getState() + r.delta_1());
        tmp = r.postProcess_T1(*p.first, *p.second);
        modified_particles.insert(tmp.begin(), tmp.end());

        p.second->setState(p.second->getState() + r.delta_2());
        tmp = r.postProcess_T2(*p.second, *p.first);
        modified_particles.insert(tmp.begin(), tmp.end());
        updated_fpl = true;
      }
    }

    // Trigger update of FixedPairList.
    if (updated_fpl)
      (*it)->fixed_pair_list_->updateParticlesStorage();
  }
}

/** Use the (A,B) list "partners" to add bonds and change the state of the
   particles accordingly.
 */
void ChemicalReaction::ApplyAR(std::set<Particle *> &modified_particles) {
  System &system = getSystemRef();
  std::set<Particle *> tmp;

  LOG4ESPP_DEBUG(theLogger, "Entering applyAR");

  // Limit number of bonds created on each of interval steps.
  longint local_bond_count = effective_pairs_.size();  // by default no limit;
  if (bond_limit_ > 0) {
    local_bond_count = 0;
    // Calculate number of bonds that could be created at this CPU.
    for (ReactionMap::const_iterator it = effective_pairs_.begin();
         it != effective_pairs_.end(); it++) {
      Particle *p1 = system.storage->lookupLocalParticle(it->first);
      Particle *p2 = system.storage->lookupLocalParticle(it->second.first);

      if (p1 && p2)
        if (!(p1->ghost() && p2->ghost()))
          local_bond_count++;
    }

    // We need to get the number of bonds on each of cpus and then redistribute
    // the correct fraction.
    std::vector<longint> global_bond_count;
    if (system.comm->rank() == 0) {
      // Collect bonds from CPUs.
      mpi::gather(*(system.comm), local_bond_count, global_bond_count, 0);

      longint total_count = std::accumulate(global_bond_count.begin(), global_bond_count.end(), 0);
      longint bonds_left = total_count > bond_limit_ ? bond_limit_ : total_count;

      std::vector<longint> global_bond_limit;
      global_bond_limit.resize(global_bond_count.size(), 0);

      // At least one bond on each of CPUs, then redistribute one by one.
      while (bonds_left > 0) {
        longint rank_index = 0;
        for (std::vector<longint>::iterator it = global_bond_count.begin();
             it != global_bond_count.end() && bonds_left > 0;
             rank_index++, ++it) {
          longint b = *it;
          if (b > 0 && global_bond_limit[rank_index] < b) {
            global_bond_limit[rank_index]++;
            bonds_left--;
          }
        }
      }

      mpi::scatter(*(system.comm), global_bond_limit, local_bond_count, 0);
    } else {
      mpi::gather(*(system.comm), local_bond_count, global_bond_count, 0);

      // Get the new value of local_bond_count.
      mpi::scatter(*(system.comm), local_bond_count, 0);
    }
  }

  for (integrator::ReactionMap::iterator it = effective_pairs_.begin();
      it != effective_pairs_.end() && local_bond_count != 0; it++) {
    boost::shared_ptr<integrator::Reaction> reaction = reaction_list_.at(it->second.second.reaction_id);

    // Change the state of A and B.
    Particle *p1 = system.storage->lookupLocalParticle(it->first);
    Particle *p2 = system.storage->lookupLocalParticle(it->second.first);

#ifdef LOG4ESPP_DEBUG_ENABLED
    if (p1 && p2) {
      LOG4ESPP_DEBUG(
          theLogger,
          "Checking pair: " << p1->id() << "(" << p1->state() << "-"
              << p2->id() << "(" << p2->state() << ") A.type="
              << p2->type() << " B.type=" << p2->type());
    }
#endif
    bool valid_state = false;

    if (p1 && p2) {
      valid_state = (reaction->isValidState_T1(*p1) && reaction->isValidState_T2(*p2));

      if (valid_state) {
        p1->setState(p1->getState() + reaction->delta_1());
        tmp = reaction->postProcess_T1(*p1, *p2);
        modified_particles.insert(p1);
        for (std::set<Particle *>::iterator pit = tmp.begin(); pit != tmp.end(); ++pit)
          modified_particles.insert(*pit);

        p2->setState(p2->getState() + reaction->delta_2());
        tmp = reaction->postProcess_T2(*p2, *p1);
        modified_particles.insert(p2);
        for (std::set<Particle *>::iterator pit = tmp.begin(); pit != tmp.end(); ++pit)
          modified_particles.insert(*pit);

        // Try to add bond.
        valid_state = reaction->fixed_pair_list_->iadd(it->first, it->second.first);
      }

      if (valid_state) {
        local_bond_count--;
      }
    }
  }

  LOG4ESPP_DEBUG(theLogger, "Leaving applyAR");
  LOG4ESPP_DEBUG(theLogger, "applyAR, modified_particles: " << modified_particles.size());
}

void ChemicalReaction::disconnect() {
  react_.disconnect();
}

void ChemicalReaction::connect() {
  react_ = integrator->aftIntV.connect(
    boost::bind(&ChemicalReaction::React, this), boost::signals2::at_front);
}

python::list ChemicalReaction::getTimers() {
  python::list ret;
  ret.append(python::make_tuple("timeComm", timeComm));
  ret.append(python::make_tuple("timeUpdateGhost", timeUpdateGhost));
  ret.append(python::make_tuple("timeApplyAR", timeApplyAR));
  ret.append(python::make_tuple("timeApplyDR", timeApplyDR));
  ret.append(python::make_tuple("timeLoopPair", timeLoopPair));

  return ret;
}

/****************************************************
** REGISTRATION WITH PYTHON
****************************************************/
void ChemicalReaction::registerPython() {
  using namespace espressopp::python;// NOLINT
  class_<ChemicalReaction, shared_ptr<ChemicalReaction>, bases<Extension> >(
      "integrator_ChemicalReaction",
      init<shared_ptr<System>, shared_ptr<VerletList>,
      shared_ptr<storage::DomainDecomposition>, shared_ptr<TopologyManager> >())
    .def("connect", &ChemicalReaction::connect)
    .def("disconnect", &ChemicalReaction::disconnect)
    .def("add_reaction", &ChemicalReaction::addReaction)
    .def("get_timers", &ChemicalReaction::getTimers)
    .add_property(
      "interval",
      &ChemicalReaction::interval,
      &ChemicalReaction::set_interval)
    .add_property(
      "nearest_mode",
      &ChemicalReaction::is_nearest,
      &ChemicalReaction::set_is_nearest)
    .add_property(
      "bond_limit",
      &ChemicalReaction::bond_limit,
      &ChemicalReaction::set_bond_limit
    );
}
}// namespace integrator
}// namespace espressopp
