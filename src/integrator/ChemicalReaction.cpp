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

#include "ChemicalReaction.hpp"

#include "storage/Storage.hpp"
#include "iterator/CellListIterator.hpp"
#include "esutil/RNG.hpp"
#include "bc/BC.hpp"
#include "storage/NodeGrid.hpp"
#include "storage/DomainDecomposition.hpp"
#include "FixDistances.hpp"

namespace espressopp {
namespace integrator {

LOG4ESPP_LOGGER(Reaction::theLogger, "Reaction");

/** Checks if the particles pair is valid. */
bool Reaction::IsValidPair(Particle &p1, Particle &p2, ParticlePair &particle_order) {
  LOG4ESPP_DEBUG(theLogger, "entering Reaction::IsValidPair, min_cutoff=" << min_cutoff_ << " cutoff_=" << cutoff_);
  if (IsValidState(p1, p2, particle_order)) {
    Real3D distance = p1.position() - p2.position();
    real distance_2 = distance.sqr();
    if (distance_2 < cutoff_sqr_ && distance_2 >= min_cutoff_sqr_ && (*rng_)() < rate_ * (*dt_) * (*interval_)) {
      LOG4ESPP_DEBUG(theLogger, "valid pair to bond " << p1.id() << "-" << p2.id() << " d2=" << distance_2);
      return true;
    }
  }
  return false;
}


/** Checks if the particles has correct state. */
bool Reaction::IsValidState(Particle &p1, Particle &p2, ParticlePair &correct_order) {
  if ((p1.res_id() == p2.res_id()) && !intramolecular_)
    return false;

  int p1_state = p1.state();
  int p2_state = p2.state();

  // States has to be always positive or zero.
  assert(p1_state >= 0 && p2_state >= 0);

  // Case when both types are the same.
  if (type_1_ == type_2_ && p1.type() == type_1_ && p1.type() == p2.type()) {
    if (p1_state >= min_state_1_ && p1_state < max_state_1_ &&
        p2_state >= min_state_2_ && p2_state < max_state_2_) {
      correct_order.first = &p1;
      correct_order.second = &p2;
      return true;
    } else if (p2_state >= min_state_1_ && p2_state < max_state_1_ &&
        p1_state >= min_state_2_ && p1_state < max_state_2_) {
      correct_order.first = &p2;
      correct_order.second = &p1;
      return true;
    }
  } else if (type_1_ != type_2_) {  // inhomogenious case.
    if ((p1.type() == type_1_) && (p2.type() == type_2_)
        && (p1_state >= min_state_1_ && p1_state < max_state_1_)
        && (p2_state >= min_state_2_ && p2_state < max_state_2_)) {
      correct_order.first = &p1;
      correct_order.second = &p2;
      return true;
    } else if ((p1.type() == type_2_) && (p2.type() == type_1_)
        && (p1_state >= min_state_2_ && p1_state < max_state_2_)
        && (p2_state >= min_state_1_ && p2_state < max_state_1_)) {
      correct_order.first = &p2;
      correct_order.second = &p1;
      return true;
    }
  }
  return false;
}

bool Reaction::IsValidStateT_1(Particle &p) {
  if (p.type() != type_1_)
    throw std::runtime_error("Particle has wrong type.");

  // States has to be always positive or zero.
  int p_state = p.state();
  assert(p_state >= 0);
  return (p_state >= min_state_1_ && p_state < max_state_1_);
}

bool Reaction::IsValidStateT_2(Particle &p) {
  if (p.type() != type_2_)
    throw std::runtime_error("Particle has wrong type.");

  // States has to be always positive or zero.
  int p_state = p.state();
  assert(p_state >= 0);
  return (p_state >= min_state_2_ && p_state < max_state_2_);
}


std::set<Particle *> Reaction::PostProcess_T1(Particle &p) {
  std::set < Particle * > output;
  std::vector < Particle * > ret;
  for (std::vector < shared_ptr < integrator::PostProcess > > ::iterator it =
                                                                             post_process_T1.begin();
      it != post_process_T1.end();
  ++it) {
    ret = (*it)->process(p);
    output.insert(ret.begin(), ret.end());
  }
  return output;
}

std::set<Particle *> Reaction::PostProcess_T2(Particle &p) {
  std::set < Particle * > output;
  std::vector < Particle * > ret;
  for (std::vector < shared_ptr < integrator::PostProcess > > ::iterator it =
                                                                             post_process_T2.begin();
      it != post_process_T2.end();
  ++it) {
    ret = (*it)->process(p);
    output.insert(ret.begin(), ret.end());
  }
  return output;
}


void Reaction::registerPython() {
  using namespace espressopp::python; //NOLINT
  class_ < Reaction, shared_ptr < integrator::Reaction > >
      ("integrator_Reaction",
          init < int, int, int, int, int, int, int, int, real, real,
          shared_ptr < FixedPairList >, bool > ())
          .add_property("type_1", &Reaction::type_1, &Reaction::set_type_1)
          .add_property("type_2", &Reaction::type_2, &Reaction::set_type_2)
          .add_property("delta_1", &Reaction::delta_1, &Reaction::set_delta_1)
          .add_property("min_state_1", &Reaction::min_state_1, &Reaction::set_min_state_1)
          .add_property("max_state_1", &Reaction::max_state_1, &Reaction::set_max_state_1)
          .add_property("delta_2", &Reaction::delta_2, &Reaction::set_delta_2)
          .add_property("min_state_2", &Reaction::min_state_2, &Reaction::set_min_state_2)
          .add_property("max_state_2", &Reaction::max_state_2, &Reaction::set_max_state_2)
          .add_property("rate", &Reaction::rate, &Reaction::set_rate)
          .add_property("cutoff", &Reaction::cutoff, &Reaction::set_cutoff)
          .add_property("min_cutoff", &Reaction::min_cutoff, &Reaction::set_min_cutoff)
          .add_property("intramolecular", &Reaction::intramolecular, &Reaction::set_intramolecular)
          .add_property("active", &Reaction::active, &Reaction::set_active)
          .def("add_postprocess", &Reaction::AddPostProcess);
}

/** DissociationReaction */
LOG4ESPP_LOGGER(DissociationReaction::theLogger, "DissociationReaction");

/** Checks if the particles pair is valid. */
bool DissociationReaction::IsValidPair(Particle &p1, Particle &p2, ParticlePair &particle_order) {
  LOG4ESPP_DEBUG(theLogger, "entering DissociationReaction::IsValidPair");
  if (IsValidState(p1, p2, particle_order)) {
    real W = (*rng_)();

    if (rate_ > 0.0) {
      Real3D distance;
      bc_->getMinimumImageVectorBox(distance, p1.position(), p2.position());
      real distance_2 = distance.sqr();

      // Break the bond when the distance exceed the cut_off with some probability.
      if (distance_2 > cutoff_sqr_ && W < rate_ * (*dt_) * (*interval_)) {
        LOG4ESPP_DEBUG(theLogger,
                       "Break the bond, " << p1.id() << "-" << p2.id()
                           << " d_2=" << distance_2 << " cutoff_sqr=" << cutoff_sqr_);
        return true;
      }
    }

    // Break the bond randomly.
    if (W < diss_rate_ * (*dt_) * (*interval_)) {
      LOG4ESPP_DEBUG(theLogger, "Break the bond randomly " << p1.id() << "-" << p2.id());
      return true;
    }
  }
  return false;
}


void DissociationReaction::registerPython() {
  LOG4ESPP_DEBUG(theLogger, "register dissociation reaction");
  using namespace espressopp::python; //NOLINT
  class_ < DissociationReaction, bases < Reaction >, shared_ptr < integrator::DissociationReaction >
      >
          ("integrator_DissociationReaction",
              init < int, int, int, int, int, int, int, int, real, real, real, shared_ptr
              < FixedPairList > > ())
              .add_property("type_1",
                            &DissociationReaction::type_1, &DissociationReaction::set_type_1)
              .add_property("type_2",
                            &DissociationReaction::type_2, &DissociationReaction::set_type_2)
              .add_property("delta_1",
                            &DissociationReaction::delta_1, &DissociationReaction::set_delta_1)
              .add_property("min_state_1",
                            &DissociationReaction::min_state_1,
                            &DissociationReaction::set_min_state_1)
              .add_property("max_state_1",
                            &DissociationReaction::max_state_1,
                            &DissociationReaction::set_max_state_1)
              .add_property("delta_2",
                            &DissociationReaction::delta_2, &DissociationReaction::set_delta_2)
              .add_property("min_state_2",
                            &DissociationReaction::min_state_2,
                            &DissociationReaction::set_min_state_2)
              .add_property("max_state_2",
                            &DissociationReaction::max_state_2,
                            &DissociationReaction::set_max_state_2)
              .add_property("rate",
                            &DissociationReaction::rate, &DissociationReaction::set_rate)
              .add_property("cutoff",
                            &DissociationReaction::cutoff, &DissociationReaction::set_cutoff)
              .add_property("diss_rate",
                            &DissociationReaction::diss_rate, &DissociationReaction::set_diss_rate)
              .add_property("active", &DissociationReaction::active, &DissociationReaction::set_active)
              .def("add_postprocess", &DissociationReaction::AddPostProcess);
}

LOG4ESPP_LOGGER(ChemicalReaction::theLogger, "ChemicalReaction");

/** ChemicalReaction part*/
ChemicalReaction::ChemicalReaction(
    shared_ptr <System> system, shared_ptr <VerletList> verletList,
    shared_ptr <storage::DomainDecomposition> domdec)
    : Extension(system),
      verlet_list_(verletList),
      domdec_(domdec) {
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
}

ChemicalReaction::~ChemicalReaction() {
  disconnect();
  LOG4ESPP_INFO(theLogger, "Destructor ChemicalReaction");
}

void ChemicalReaction::Initialize() {
  LOG4ESPP_INFO(theLogger, "Init ChemicalReaction");
}

/** Adds the chemical reaction to the list of reactions */
void ChemicalReaction::AddReaction(boost::shared_ptr <integrator::Reaction> reaction) {
  reaction->set_dt(dt_);
  reaction->set_interval(interval_);
  reaction->set_rng(rng_);
  bc::BC &bc = *getSystemRef().bc;
  reaction->set_bc(&bc);

  if (!reaction->reverse()) {
    // If VL cutoff is smaller than reaction, increase it.
    if (reaction->cutoff() > current_cutoff_) {
      LOG4ESPP_INFO(theLogger, "VL cutoff is extended to match with reaction cutoff");
      verlet_list_->setVerletCutoff(reaction->cutoff());
    }
    LOG4ESPP_INFO(theLogger, "Add reaction, rate=" << reaction->rate() * (*dt_) * (*interval_));
    reaction_list_.push_back(reaction);
  } else {
    LOG4ESPP_INFO(theLogger, "Add reverse reaction, rate=" << reaction->rate() * (*dt_) * (*interval_));
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

  // Copy particle tuple lists of reverse reactions.
  // CopyPairList();

  *dt_ = integrator->getTimeStep();

  potential_pairs_.clear();
  effective_pairs_.clear();
  // loop over VL pairs
  for (PairList::Iterator it(verlet_list_->getPairs()); it.isValid(); ++it) {
    Particle &p1 = *it->first;
    Particle &p2 = *it->second;
    int reaction_idx_ = 0;

    for (ReactionList::iterator it = reaction_list_.begin();
         it != reaction_list_.end(); ++it, ++reaction_idx_) {
      if (!(*it)->active())
        continue;
      ParticlePair p;
      if ((*it)->IsValidPair(p1, p2, p)) {
        potential_pairs_.insert(
            std::make_pair(p.first->id(), std::make_pair(p.second->id(), reaction_idx_)));
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
  std::set < Particle * > modified_particles;
  ApplyDR(modified_particles);
  ApplyAR(modified_particles);
//  UpdateGhost(modified_particles);
//  modified_particles.clear();
  // Update the ghost particles.
  UpdateGhost(modified_particles);
  LOG4ESPP_DEBUG(theLogger, "Finished react()");
  LOG4ESPP_DEBUG(theLogger, "Leaving react()");
}

/** Performs two-way parallel communication to consolidate mm between
 neighbours. The parallel scheme is taken from
 storage::DomainDecomposition::doGhostCommunication
 */
void ChemicalReaction::SendMultiMap(integrator::ReactionMap &mm) {  //NOLINT
  LOG4ESPP_DEBUG(theLogger, "Entering sendMultiMap");

  InBuffer in_buffer_0(*getSystem()->comm);
  InBuffer in_buffer_1(*getSystem()->comm);
  OutBuffer out_buffer(*getSystem()->comm);
  const storage::NodeGrid &node_grid = domdec_->getNodeGrid();

  // Prepare out buffer with the reactions that potential will happen on this node.
  out_buffer.reset();
  in_buffer_0.reset();
  in_buffer_1.reset();

  // Fill out_buffer from mm.
  int tmp = mm.size();
  int a, b, c;
  out_buffer.write(tmp);
  for (integrator::ReactionMap::iterator it = mm.begin(); it != mm.end();
       it++) {
    a = it->first;  // particle id
    b = it->second.first;  // particle id
    c = it->second.second;  // reaction id
    out_buffer.write(a);
    out_buffer.write(b);
    out_buffer.write(c);
  }

  /* direction loop: x, y, z.
   Here we could in principle build in a one sided ghost
   communication, simply by taking the lr loop only over one
   value. */

  int data_length, idx_a, idx_b, reaction_idx, direction_size;
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
        } else {
          in_buffer_1.read(idx_a);
          in_buffer_1.read(idx_b);
          in_buffer_1.read(reaction_idx);
        }
        mm.insert(std::make_pair(idx_a, std::make_pair(idx_b, reaction_idx)));
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
void ChemicalReaction::UpdateGhost(const std::set<Particle *> &modified_particles) {  //NOLINT
  LOG4ESPP_DEBUG(theLogger, "Entering UpdateGhost");

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
  longint p_id, p_type, p_res_id;
  real p_mass, p_q;
  out_buffer.write(data_length);
  for (std::set<Particle *>::const_iterator it = modified_particles.begin();
       it != modified_particles.end();
       ++it) {
    p_id = (*it)->id();
    p_type = (*it)->type();
    p_mass = (*it)->mass();
    p_q = (*it)->q();
    p_res_id = (*it)->res_id();
    out_buffer.write(p_id);
    out_buffer.write(p_type);
    out_buffer.write(p_mass);
    out_buffer.write(p_q);
    out_buffer.write(p_res_id);
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
        } else {
          in_buffer_1.read(p_id);
          in_buffer_1.read(p_type);
          in_buffer_1.read(p_mass);
          in_buffer_1.read(p_q);
          in_buffer_1.read(p_res_id);
        }
        // Update the ghost particle data on neighbour CPUs.
        particle = system.storage->lookupLocalParticle(p_id);
        if (particle != NULL) {
          LOG4ESPP_DEBUG(theLogger, "Update particle data");
          particle->setType(p_type);
          particle->setMass(p_mass);
          particle->setQ(p_q);
          particle->setResId(p_res_id);
        }
      }
    }
    LOG4ESPP_DEBUG(theLogger, "Leaving unpack");
  }
  LOG4ESPP_DEBUG(theLogger, "Leaving UpdateGhost");
}

/** Given a multimap mm with several pairs (id1,id2), keep only one pair for
 each id1 and return it in place. In addition, only pairs for which
 id1 is local are kept.
 */
void ChemicalReaction::UniqueA(integrator::ReactionMap &potential_candidates) {  //NOLINT
  LOG4ESPP_DEBUG(theLogger, "UniqueA");
  System &system = getSystemRef();
  integrator::ReactionMap unique_list_of_candidates;
  boost::unordered_set <longint> a_indexes;

  Particle *p = NULL;
  unique_list_of_candidates.clear();
  // Gets the list of indexes of particle a. Gets only real particles, skip ghost.
  for (integrator::ReactionMap::iterator it = potential_candidates.begin();
       it != potential_candidates.end(); ++it) {
    p = system.storage->lookupLocalParticle(it->first);
    if (p == NULL || p->ghost())
      continue;
    a_indexes.insert(it->first);
  }

  // For each active idx1, pick a partner
  if (a_indexes.size() > 0) {
    int idx_a;
    real max_reaction_rate;
    // rate => idx_b, reaction_id
    boost::unordered_multimap <real, std::pair<longint, int> > rate_idx_b;
    boost::unordered_multimap < real, std::pair < longint, int > > ::local_iterator
    idx_b_reaction_id;

    // Iterators for the equal_range.
    std::pair <integrator::ReactionMap::iterator,
    integrator::ReactionMap::iterator> candidates_b;
    // Iterate over the ids of particle A, looking for the particle B
    for (boost::unordered_set<longint>::iterator it = a_indexes.begin();
         it != a_indexes.end(); ++it) {
      idx_a = *it;

      // Select all possible candidates
      candidates_b = potential_candidates.equal_range(idx_a);

      // Group the candidates by the reaction rate.
      max_reaction_rate = -1;
      rate_idx_b.clear();
      for (integrator::ReactionMap::iterator jt = candidates_b.first;
           jt != candidates_b.second; ++jt) {
        boost::shared_ptr <integrator::Reaction> reaction = reaction_list_.at(jt->second.second);
        real reaction_rate = reaction->rate();
        if (reaction_rate > max_reaction_rate) {
          max_reaction_rate = reaction_rate;
        }
        // rate => (idx_b, reaction_id)
        rate_idx_b.insert(
            std::make_pair(
                reaction_rate,
                std::make_pair(jt->second.first, jt->second.second)));
      }

      // Found reaction with the maximum rate. If there are several candidates with the same
      // rate, then we choose randomly.
      if (max_reaction_rate != -1) {
        int bucket_size = rate_idx_b.count(max_reaction_rate);
        // Pick up random number in given range.
        int pick_offset = (*rng_)(bucket_size);
        idx_b_reaction_id = rate_idx_b.begin(
            rate_idx_b.bucket(max_reaction_rate));
        std::advance(idx_b_reaction_id, pick_offset);

        unique_list_of_candidates.insert(
            std::make_pair(idx_a, idx_b_reaction_id->second));
      }
    }
  }
  potential_candidates = unique_list_of_candidates;
}

/** Given a multimap mm with several pairs (id1,id2), keep only one pair for
 each id2 and return it in place. In addition, only pairs for which
 id2 is local are kept.
 */
void ChemicalReaction::UniqueB(integrator::ReactionMap &potential_candidates,  //NOLINT
                               integrator::ReactionMap &effective_candidates) {  //NOLINT
  LOG4ESPP_DEBUG(theLogger, "UniqueB");

  typedef boost::unordered_set <longint> Indexes;
  typedef boost::unordered_multimap <real, std::pair<longint, int> > RateParticleIdx;

  System &system = getSystemRef();
  Indexes b_indexes;
  integrator::ReactionMap reverse_candidates;

  Particle *p = NULL;
  effective_candidates.clear();

  // Collect the b particle pairs. REQOPT
  for (integrator::ReactionMap::iterator it = potential_candidates.begin();
       it != potential_candidates.end(); ++it) {
    p = system.storage->lookupLocalParticle(it->second.first);
    if (p == NULL || p->ghost())
      continue;

    b_indexes.insert(it->second.first);
    reverse_candidates.insert(
        std::make_pair(it->second.first,
                       std::make_pair(it->first, it->second.second)));
  }

  if (b_indexes.size() > 0) {
    int idx_b;
    real max_reaction_rate;
    // rate => idx_a, reaction_id
    RateParticleIdx rate_idx_a;
    RateParticleIdx::local_iterator idx_a_reaction_id;
    std::pair <integrator::ReactionMap::iterator,
    integrator::ReactionMap::iterator> candidates_a;

    for (Indexes::iterator it = b_indexes.begin(); it != b_indexes.end(); ++it) {
      idx_b = *it;

      candidates_a = reverse_candidates.equal_range(idx_b);
      max_reaction_rate = -1;
      rate_idx_a.clear();
      for (integrator::ReactionMap::iterator jt = candidates_a.first;
           jt != candidates_a.second; ++jt) {
        boost::shared_ptr <integrator::Reaction> reaction = reaction_list_.at(jt->second.second);
        real reaction_rate = reaction->rate();
        if (reaction_rate > max_reaction_rate) {
          max_reaction_rate = reaction_rate;
        }
        rate_idx_a.insert(
            std::make_pair(
                reaction_rate,
                std::make_pair(jt->second.first, jt->second.second)));
      }
      // Found reaction with the maximum rate. If there are several candidates
      // then select randomly.
      if (max_reaction_rate > -1) {
        int bucket_size = rate_idx_a.count(max_reaction_rate);
        int pick_offset = (*rng_)(bucket_size);
        idx_a_reaction_id = rate_idx_a.begin(
            rate_idx_a.bucket(max_reaction_rate));
        std::advance(idx_a_reaction_id, pick_offset);
        effective_candidates.insert(
            std::make_pair(
                idx_a_reaction_id->second.first,
                std::make_pair(idx_b, idx_a_reaction_id->second.second)));
      }
    }
  }
}


/** ApplyDR to remove bonds and change the state of the particles
 * accordingly
 */
void ChemicalReaction::ApplyDR(std::set < Particle * > &modified_particles) {
  LOG4ESPP_DEBUG(theLogger, "Entering applyDR");

  // Iterate over reverse reaction. For every reaction, iterate overy particles pairs and
  // decide to remove or keep bond.
  std::set < Particle * > tmp;
  for (ReactionList::iterator it = reverse_reaction_list_.begin();
       it != reverse_reaction_list_.end(); ++it) {
    integrator::Reaction &r = **it;
    if (!r.active())
      continue;
    bool updated_fpl = false;
    for (FixedPairList::PairList::Iterator itp(*((*it)->fixed_pair_list_)); itp.isValid(); ++itp) {
      Particle &p1 = *itp->first;
      Particle &p2 = *itp->second;
      ParticlePair p;
      if ((*it)->IsValidPair(p1, p2, p)) {
        // Ok we will remove this pair.
        if (!p1.ghost() && !p1.ghost())
          (*it)->fixed_pair_list_->remove(p1.id(), p2.id());
        else if (p2.ghost())
          (*it)->fixed_pair_list_->remove(p1.id(), p2.id());
        else
          (*it)->fixed_pair_list_->remove(p2.id(), p1.id());
        // Run postprocess.
        p.first->setState(p.first->getState() + r.delta_1());
        tmp = r.PostProcess_T1(*p.first);
        modified_particles.insert(tmp.begin(), tmp.end());

        p.second->setState(p.second->getState() + r.delta_2());
        tmp = r.PostProcess_T2(*p.second);
        modified_particles.insert(tmp.begin(), tmp.end());
        updated_fpl = true;
      }
    }
    // Trigger update of FixedPairList.
    if (updated_fpl)
      (*it)->fixed_pair_list_->onParticlesChanged();
  }
}


/** Use the (A,B) list "partners" to add bonds and change the state of the
 particles accordingly.
 */
void ChemicalReaction::ApplyAR(std::set < Particle * > &modified_particles) {
  System &system = getSystemRef();

  std::set < Particle * > tmp;

  LOG4ESPP_DEBUG(theLogger, "Entering applyAR");

  for (integrator::ReactionMap::iterator it = effective_pairs_.begin();
       it != effective_pairs_.end(); it++) {
    boost::shared_ptr <integrator::Reaction> reaction = reaction_list_.at(it->second.second);
    // Change the state of A and B.
    Particle *p1 = system.storage->lookupLocalParticle(it->first);
    Particle *p2 = system.storage->lookupLocalParticle(it->second.first);
    LOG4ESPP_DEBUG(
        theLogger,
        "Checking pair: " << p1->id() << "(" << p1->state() << "-"
            << p2->id() << "(" << p2->state() << ") A.type="
            << p2->type() << " B.type=" << p2->type());

    bool valid_state = true;
    if (p1 != NULL) {
      if (reaction->IsValidStateT_1(*p1)) {
        p1->setState(p1->getState() + reaction->delta_1());
        tmp = reaction->PostProcess_T1(*p1);
        for (std::set<Particle *>::iterator pit = tmp.begin(); pit != tmp.end(); ++pit)
          modified_particles.insert(*pit);
      } else {
        valid_state = false;
      }
    }
    if (p2 != NULL) {
      if (reaction->IsValidStateT_2(*p2)) {
        p2->setState(p2->getState() + reaction->delta_2());
        tmp = reaction->PostProcess_T2(*p2);
        for (std::set<Particle *>::iterator pit = tmp.begin(); pit != tmp.end(); ++pit)
          modified_particles.insert(*pit);
      } else {
        valid_state = false;
      }
    }
    /** Make sense only if both particles exists here, otherwise waste of CPU time. */
    if (p1 != NULL && p2 != NULL && valid_state) {
      LOG4ESPP_DEBUG(theLogger, "adding pair " << it->first << "-" << it->second.first);
      reaction->fixed_pair_list_->iadd(it->first, it->second.first);
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


/****************************************************
 ** REGISTRATION WITH PYTHON
 ****************************************************/
void ChemicalReaction::registerPython() {
  using namespace espressopp::python; //NOLINT
  class_ < ChemicalReaction, shared_ptr < ChemicalReaction >, bases < Extension > > (
      "integrator_ChemicalReaction",
          init < shared_ptr < System >, shared_ptr < VerletList >,
          shared_ptr < storage::DomainDecomposition > > ())
      .def("connect", &ChemicalReaction::connect)
      .def("disconnect", &ChemicalReaction::disconnect)
      .def("add_reaction", &ChemicalReaction::AddReaction)
      .add_property(
          "interval",
          &ChemicalReaction::interval,
          &ChemicalReaction::set_interval);
}

LOG4ESPP_LOGGER(PostProcess::theLogger, "PostProcess");

/** PostProcess methods */
void PostProcess::registerPython() {
  using namespace espressopp::python;  //NOLINT
  class_ < PostProcess, shared_ptr < integrator::PostProcess >, boost::noncopyable >
      ("integrator_PostProcess", no_init);
}

/** Adds new change property definition. */
void PostProcessChangeProperty::AddChangeProperty(
    int type_id,
    boost::shared_ptr <ParticleProperties> new_property) {
  std::pair<TypeParticlePropertiesMap::iterator, bool> ret;
  ret = type_properties_.insert(
      std::pair < int, boost::shared_ptr < ParticleProperties > > (type_id, new_property));
  if (ret.second == false)
    throw std::runtime_error("Requested type already exists. To replace please remove it firstly");
}

/** Removes change property definition. */
void PostProcessChangeProperty::RemoveChangeProperty(int type_id) {
  int remove_elements = type_properties_.erase(type_id);
  if (remove_elements == 0) {
    throw std::runtime_error("Invalid type.");
  }
}

/** Post process after pairs were added.
 *
 * In this case method will update the properties of the particles.
 * */
std::vector<Particle *> PostProcessChangeProperty::process(Particle &p1) {
  TypeParticlePropertiesMap::iterator it;
  std::vector < Particle *> mod_particles;
  LOG4ESPP_DEBUG(theLogger, "Entering PostProcessChangeProperty::process()");
  // Process particle p1.
  it = type_properties_.find(p1.type());
  bool mod = false;
  LOG4ESPP_DEBUG(theLogger, "type " << it->second->type);
  if (it != type_properties_.end()) {
    if (it->second->type != NULL) {
      p1.setType(it->second->type);
      mod = true;
    }
    if (it->second->mass != NULL) {
      p1.setMass(it->second->mass);
      mod = true;
    }
    if (it->second->q != NULL) {
      p1.setQ(it->second->q);
      mod = true;
    }
    LOG4ESPP_DEBUG(theLogger, "Modified particle A: " << p1.id());
  }
  if (mod)
    mod_particles.push_back(&(p1));
  return mod_particles;
}

void PostProcessChangeProperty::registerPython() {
  using namespace espressopp::python;  //NOLINT

  class_ < PostProcessChangeProperty, bases < integrator::PostProcess >,
      boost::shared_ptr < integrator::PostProcessChangeProperty > >
          ("integrator_PostProcessChangeProperty", init<>())
              .def("add_change_property", &PostProcessChangeProperty::AddChangeProperty)
              .def("remove_change_property", &PostProcessChangeProperty::RemoveChangeProperty);
}

}  // namespace integrator
}  // namespace espressopp
