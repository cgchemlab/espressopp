/*
  Copyright (C) 2016
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

#include "python.hpp"
#include "ATRPActivator.hpp"

#include "esutil/RNG.hpp"
#include "boost/range/algorithm.hpp"

namespace espressopp {
namespace integrator {

/** ATRPActivator */
LOG4ESPP_LOGGER(ATRPActivator::theLogger, "ATRPActivator");

ATRPActivator::ATRPActivator(
    shared_ptr<System> system, longint interval, longint num_per_interval, real ratio_activator,
    real ratio_deactivator, real delta_catalyst, real k_activate, real k_deactivate)
    : Extension(system), interval_(interval), num_per_interval_(num_per_interval),
      ratio_activator_(ratio_activator), ratio_deactivator_(ratio_deactivator), delta_catalyst_(delta_catalyst),
      k_activate_(k_activate), k_deactivate_(k_deactivate){
  LOG4ESPP_INFO(theLogger, "ATRPActivator constructed");

  // Set RNG.
  rng_ = system->rng;
}

void ATRPActivator::disconnect() {
  sig_aftIntV.disconnect();
}

void ATRPActivator::connect() {
  sig_aftIntV = integrator->aftIntV.connect(boost::bind(&ATRPActivator::updateParticles, this));
}

void ATRPActivator::addReactiveCenter(longint type_id,
                                      longint min_state,
                                      longint max_state,
                                      shared_ptr<ParticleProperties> pp,
                                      longint delta_state) {
  if (min_state > max_state)
    throw std::runtime_error("Min_state > max_state");

  if (species_map_.count(type_id) != 0) {
    std::pair<SpeciesMap::iterator, SpeciesMap::iterator> equalRange = species_map_.equal_range(type_id);
    std::vector<std::pair<longint, longint> > ranges;
    ranges.push_back(std::make_pair(min_state, max_state));  // put the possible range
    for (SpeciesMap::iterator it = equalRange.first; it != equalRange.second; ++it) {
      ranges.push_back(std::make_pair(it->second.min_state, it->second.max_state));
    }
    std::sort(ranges.begin(), ranges.end());  // sort by first element, min_state
    bool overlaps = false;
    for (int i = 0; i < ranges.size() - 1; i++) {
      if (ranges[i].second > ranges[i+1].first) {
        overlaps = true;
      }
    }
    if (overlaps)
      throw std::runtime_error("Min/max state overlaps");
  }
  species_map_.insert(std::make_pair(type_id, ReactiveCenter(min_state, max_state, delta_state, pp)));
}

void ATRPActivator::updateParticles() {
  if (integrator->getStep() % (interval_) != 0)
    return;

  // Collect on every process, list of pids of given types, active or dormant species.
  std::vector<longint> local_type_pids;
  System &system = getSystemRef();
  CellList cl = system.storage->getRealCells();

  for (espressopp::iterator::CellListIterator cit(cl); !cit.isDone(); ++cit) {
    Particle &p = *cit;
    bool found = false;
    if (species_map_.count(p.type()) != 0) {
      // Check if the particle is in the given state.
      std::pair<SpeciesMap::iterator, SpeciesMap::iterator> equalRange = species_map_.equal_range(p.type());
      for (SpeciesMap::iterator it = equalRange.first; it != equalRange.second && !found; ++it) {
        if (p.state() >= it->second.min_state && p.state() < it->second.max_state) {
          local_type_pids.push_back(p.id());
          found = true;
        }
      }
    }
  }

  // Every CPU sends list of particle ids to root CPU.
  std::vector<std::vector<longint> > global_type_pids;
  std::vector<longint> selected_pids;

  if (system.comm->rank() == 0) {
    mpi::gather(*(system.comm), local_type_pids, global_type_pids, 0);

    // Root CPU select randomly :num_per_interval particles.
    std::vector<longint> all_pids;
    for (std::vector<std::vector<longint> >::iterator it = global_type_pids.begin();
         it != global_type_pids.end(); ++it) {
      all_pids.insert(all_pids.end(), it->begin(), it->end());
    }
    shared_ptr<esutil::RNG> rng = system.rng;
    boost::range::random_shuffle(all_pids, *(rng));
    for (int n = 0; n < num_per_interval_ && n < all_pids.size(); n++) {
      selected_pids.push_back(all_pids[n]);
    }
  } else {
    mpi::gather(*(system.comm), local_type_pids, global_type_pids, 0);
  }

  // Broadcast selected particles id to all CPUs.
  // TODO(jakub): instead of broadcast all selected pids, send only to CPUs that posses those particles.
  mpi::broadcast(*(system.comm), selected_pids, 0);

  // On every CPU, iterate over list of pids. If particle is on this CPU then process it.
  std::vector<Particle *> modified_particles;
  for (std::vector<longint>::iterator it = selected_pids.begin(); it != selected_pids.end(); ++it) {
    Particle *p = system.storage->lookupRealParticle(*it);  // check if particle is on this CPU.
    if (p) {
      bool found = false;
      std::pair<SpeciesMap::iterator, SpeciesMap::iterator> equalRange = species_map_.equal_range(p->type());
      for (SpeciesMap::iterator it = equalRange.first; it != equalRange.second && !found; ++it) {
        if (p->state() == it->second.min_state) {
          real W = (*rng_)();
          if (W < ratio_deactivator_*k_deactivate_) {
            p->state() += it->second.delta_state;  // update chemical state.
            it->second.new_property->updateParticleProperties(p); // update particle property.
            modified_particles.push_back(p);
            ratio_deactivator_-=delta_catalyst_;
            ratio_activator_+=delta_catalyst_;
          }
          found = true;
        }
        else if (p->state() > it->second.min_state && p->state() < it->second.max_state) {
          real W = (*rng_)();
          if (W < ratio_activator_*k_activate_) {
            p->state() += it->second.delta_state;  // update chemical state.
            it->second.new_property->updateParticleProperties(p); // update particle property.
            modified_particles.push_back(p);
            ratio_activator_-=delta_catalyst_;
            ratio_deactivator_+=delta_catalyst_;
          }
          found = true;
        }
      }
      if (!found) {
        std::stringstream ss;
        ss << "Wrong particle on rank: " << system.comm->rank()
           << " particle in state: " << p->state()
           << " id: " << p->id()
           << " type: " << p->type();
        throw std::runtime_error(ss.str());
      }

    }
  }
  // Synchronize all processess
  (*system.comm).barrier();

  // Update neighbour ghosts
  updateGhost(modified_particles);
}


/** Performs two-way parallel communication to update the ghost particles.
 * The parallel scheme is taken from
 * storage::DomainDecomposition::doGhostCommunication
 */
void ATRPActivator::updateGhost(const std::vector<Particle *> &modified_particles) {// NOLINT
  LOG4ESPP_DEBUG(theLogger, "Entering updateGhost");

  int kCrCommTag = 0x6b;

  System &system = getSystemRef();

  InBuffer in_buffer_0(*getSystem()->comm);
  InBuffer in_buffer_1(*getSystem()->comm);
  OutBuffer out_buffer(*getSystem()->comm);
  const storage::NodeGrid
      &node_grid = boost::static_pointer_cast<storage::DomainDecomposition>(system.storage)->getNodeGrid();

  // Prepare out buffer with the reactions that potential will happen on this node.
  out_buffer.reset();
  in_buffer_0.reset();
  in_buffer_1.reset();

  // Fill out_buffer from the particles properties.
  longint data_length = modified_particles.size();
  longint p_id, p_type, p_state;

  out_buffer.write(data_length);

  for (std::vector<Particle *>::const_iterator it = modified_particles.begin(); it != modified_particles.end(); ++it) {
    p_id = (*it)->id();
    p_type = (*it)->type();
    p_state = (*it)->state();
    out_buffer.write(p_id);
    out_buffer.write(p_type);
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
          in_buffer_0.read(p_state);
        } else {
          in_buffer_1.read(p_id);
          in_buffer_1.read(p_type);
          in_buffer_1.read(p_state);
        }

        // Update the ghost particle data on neighbour CPUs.
        particle = system.storage->lookupLocalParticle(p_id);

        if (particle) {
          LOG4ESPP_DEBUG(theLogger, "Update particle data");
          particle->setType(p_type);
          particle->setState(p_state);
        }
      }
    }

    LOG4ESPP_DEBUG(theLogger, "Leaving unpack");
  }

  LOG4ESPP_DEBUG(theLogger, "Leaving updateGhost");
}

void ATRPActivator::registerPython() {
  using namespace espressopp::python;  // NOLINT

  class_<ATRPActivator, shared_ptr<ATRPActivator>, bases<Extension> >
      ("integrator_ATRPActivator", init<shared_ptr<System>, longint, longint, real, real, real, real, real>())
      .def("add_reactive_center", &ATRPActivator::addReactiveCenter)
      .def("update_particles", &ATRPActivator::updateParticles)
      .def("connect", &ATRPActivator::connect)
      .def("disconnect", &ATRPActivator::disconnect);
}


}  // end namespace integrator
}  // end namespace espressopp
