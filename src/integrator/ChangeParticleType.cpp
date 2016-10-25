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
#include "ChangeParticleType.hpp"

#include "esutil/RNG.hpp"
#include "boost/range/algorithm.hpp"

namespace espressopp {
namespace integrator {

/** ChangeParticleType */
LOG4ESPP_LOGGER(ChangeParticleType::theLogger, "ChangeParticleType");

ChangeParticleType::ChangeParticleType(
    shared_ptr<System> system, longint interval, longint num_per_interval, longint old_type, longint new_type)
    : Extension(system), interval_(interval), num_per_interval_(num_per_interval),
      old_type_(old_type), new_type_(new_type) {
  LOG4ESPP_INFO(theLogger, "ChangeParticleType constructed");
}

void ChangeParticleType::disconnect() {
  sig_aftIntV.disconnect();
}

void ChangeParticleType::connect() {
  sig_aftIntV = integrator->aftIntV.connect(boost::bind(&ChangeParticleType::updateParticles, this));
}

void ChangeParticleType::updateParticles() {
  if (integrator->getStep() % (interval_) != 0)
    return;

  // Collect on every process, list of pids of given type.
  std::vector<longint> local_type_pids;
  System &system = getSystemRef();
  CellList cl = system.storage->getRealCells();
  for (espressopp::iterator::CellListIterator cit(cl); !cit.isDone(); ++cit) {
    Particle &p = *cit;
    if (p.type() == old_type_) {
      local_type_pids.push_back(p.id());
    }
  }

  // Communicate list to master node (rank = 0)
  std::vector<std::vector<longint> > global_type_pids;
  std::vector<longint> selected_pids;

  if (system.comm->rank() == 0) {
    mpi::gather(*(system.comm), local_type_pids, global_type_pids, 0);

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

  // Broadcast pids
  mpi::broadcast(*(system.comm), selected_pids, 0);

  std::set<Particle *> modified_particles;
  for (std::vector<longint>::iterator it = selected_pids.begin(); it != selected_pids.end(); ++it) {
    Particle *p = system.storage->lookupRealParticle(*it);
    if (p) {
      if (p->type() == old_type_) {
        p->type() = new_type_;
        modified_particles.insert(p);
      }
    }
  }
  // Synchronize all processess
  (*system.comm).barrier();

  // Update neighbour ghost
  updateGhost(modified_particles);
}

/** Performs two-way parallel communication to update the ghost particles.
 * The parallel scheme is taken from
 * storage::DomainDecomposition::doGhostCommunication
 */
void ChangeParticleType::updateGhost(const std::set<Particle *> &modified_particles) {// NOLINT
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
  longint p_id, p_type;

  out_buffer.write(data_length);

  for (std::set<Particle *>::const_iterator it = modified_particles.begin(); it != modified_particles.end(); ++it) {
    p_id = (*it)->id();
    p_type = (*it)->type();
    out_buffer.write(p_id);
    out_buffer.write(p_type);
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
        } else {
          in_buffer_1.read(p_id);
          in_buffer_1.read(p_type);
        }

        // Update the ghost particle data on neighbour CPUs.
        particle = system.storage->lookupLocalParticle(p_id);

        if (particle) {
          LOG4ESPP_DEBUG(theLogger, "Update particle data");
          particle->setType(p_type);
        }
      }
    }

    LOG4ESPP_DEBUG(theLogger, "Leaving unpack");
  }

  LOG4ESPP_DEBUG(theLogger, "Leaving updateGhost");
}

void ChangeParticleType::registerPython() {
  using namespace espressopp::python;  // NOLINT

  class_<ChangeParticleType, shared_ptr<ChangeParticleType>, bases<Extension> >
      ("integrator_ChangeParticleType", init<shared_ptr<System>,
                                             longint,
                                             longint,
                                             longint,
                                             longint>())
      .def("update_particles", &ChangeParticleType::updateParticles)
      .def("connect", &ChangeParticleType::connect)
      .def("disconnect", &ChangeParticleType::disconnect);
}

}  // end namespace integrator
}  // end namespace espressopp
