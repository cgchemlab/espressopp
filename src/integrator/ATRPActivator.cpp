/*
  Copyright (C) 2016
      Jakub Krajniak (jkrajniak at gmail.com)
      Zidan Zhang (zidan.zhang at kuleuven.be)

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

#include <fstream>
#include "python.hpp"
#include "ATRPActivator.hpp"

#include "esutil/RNG.hpp"
//#include "boost/range/algorithm.hpp"
#include "boost/filesystem/operations.hpp"

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
  stats_filename_ = "atrp_stats.dat";

  extensionOrder = Extension::beforeReaction;
  max_property_id_ = 0;

  select_from_all_ = true;

  resetTimers();
}

void ATRPActivator::disconnect() {
  sig_aftIntV.disconnect();
}

void ATRPActivator::connect() {
  sig_aftIntV = integrator->aftIntV.connect(extensionOrder, boost::bind(&ATRPActivator::updateParticles, this));
}

void ATRPActivator::addReactiveCenter(longint type_id,
                                      longint state,
                                      bool is_activator,
                                      shared_ptr<TopologyParticleProperties> pp,
                                      longint delta_state) {
  LOG4ESPP_DEBUG(theLogger, "ATRPActivator::addReactiveCenter type_id: " << type_id << " state:" << state
                            << " is_activator:" << is_activator << " delta_state:" << delta_state);
  if (species_map_.find(std::make_pair(type_id, state)) != species_map_.end()) {
    std::cout << "type_id=" << type_id << " state=" << state << " already defined!";
    throw std::runtime_error("ATRPActivator: duplicated entry in species_map_");
  }
  species_map_.insert(std::make_pair(
      std::make_pair(type_id, state),
      ReactiveCenter(state, is_activator, delta_state, max_property_id_)));
  property_map_.insert(std::make_pair(max_property_id_, pp));
  max_property_id_++;
}

void ATRPActivator::updateParticles() {
  real time0 = wallTimer.getElapsedTime();
  if (integrator->getStep() % (interval_) != 0)
    return;

  LOG4ESPP_INFO(theLogger, "ATRPActivator::updateParticles begin");

  System &system = getSystemRef();

  stats_k_activator.push_back(integrator->getStep());
  stats_k_activator.push_back(ratio_activator_);
  stats_k_activator.push_back(ratio_deactivator_);

  if (system.comm->rank() == 0) {
    if (!boost::filesystem::exists(stats_filename_.c_str())) {
      std::fstream fs(stats_filename_.c_str(), std::fstream::out);
      fs << "# step r_activator r_deactivator" << std::endl;
    }
    std::fstream fs(stats_filename_.c_str(), std::fstream::out | std::fstream::app);
    fs << integrator->getStep() << " "
       << std::scientific << ratio_activator_ << " "
       << std::scientific << ratio_deactivator_ << std::endl;
    fs.close();
  }

  // Collect on every process, list of pids of given types, active or dormant species.
  std::vector<longint> local_type_pids_state;
  CellList cl = system.storage->getRealCells();

  longint my_N = 0;

  // TODO(jakub): this is expensive part
  for (espressopp::iterator::CellListIterator cit(cl); !cit.isDone(); ++cit) {
    Particle &p = *cit;
    bool found = false;
    std::pair<longint, longint> key(p.type(), p.state());
    if (species_map_.count(key) != 0) {
      local_type_pids_state.push_back(p.id());
      local_type_pids_state.push_back(p.type());
      local_type_pids_state.push_back(p.state());
    }
    my_N++;
  }

  // Every CPU sends list of particle ids to root CPU.
  std::vector<std::vector<longint> > global_type_pids;
  std::vector<longint> selected_pids;
  longint total_N = 0;

  longint my_rank = system.comm->rank();

  if (my_rank == 0) {
    // Get the total number of particles.
    mpi::reduce(*(system.comm), my_N, total_N, std::plus<longint>(), 0);
    mpi::gather(*(system.comm), local_type_pids_state, global_type_pids, 0);

    // Root CPU select randomly :num_per_interval particles. Flat the list.
    std::map<longint, ATRPParticleP*> all_pids_state;
    for (std::vector<std::vector<longint > >::iterator it = global_type_pids.begin();
            it != global_type_pids.end(); ++it) {
      for (std::vector<longint>::iterator itt = it->begin(); itt != it->end();) {
        longint p_id = *(itt++);
        longint p_type = *(itt++);
        longint p_state = *(itt++);
        ATRPParticleP *atrpParticleP = new ATRPParticleP(p_id, p_type, p_state);
        all_pids_state.insert(std::make_pair(p_id, atrpParticleP));
      }
    }

    longint num_particles = all_pids_state.size();
    // MC step
    longint accept = 0;
    longint act = 0;
    longint dact = 0;

    // If the switch is false then we limited this loop
    if (!select_from_all_)
      total_N = num_particles;
    std::cout << "select_from_all_: " << select_from_all_ << std::endl;

    longint idx = 0;
    ATRPParticleP *atrpParticleP;
    std::map<longint, ATRPParticleP*>::iterator it_all_pids_state;

    for (int n = 0; n < total_N; n++) {  // internal MC trial
      if (select_from_all_) {
        idx = (*rng_)(total_N + 1);
        if (all_pids_state.find(idx) == all_pids_state.end())
          continue;

        accept++;
        atrpParticleP = all_pids_state.at(idx);
      } else {
        idx = (*rng_)(total_N);
        it_all_pids_state = all_pids_state.begin();
        std::advance(it_all_pids_state, idx);
        atrpParticleP = it_all_pids_state->second;
        accept++;
      }

      if (atrpParticleP == NULL || !atrpParticleP->init)
        throw std::runtime_error("wrong atrpParticleP");

      longint p_id = atrpParticleP->p_id;
      longint p_type = atrpParticleP->p_type;
      longint p_state = atrpParticleP->p_state;

      if (species_map_.count(std::make_pair(p_type, p_state)) == 1) {
        ReactiveCenter rc = species_map_.at(std::make_pair(p_type, p_state));
        real W = (*rng_)();
        shared_ptr<TopologyParticleProperties> prop = property_map_.at(rc.property_id);
        if (rc.is_activator) {
          if (W < ratio_deactivator_*k_deactivate_) {
            ratio_deactivator_ -= delta_catalyst_;
            ratio_activator_ += delta_catalyst_;
            atrpParticleP->updated = true;
            atrpParticleP->property_id = rc.property_id;
            atrpParticleP->p_state += rc.delta_state;
            atrpParticleP->p_type = prop->type();
            dact++;
          }
        } else {
          if (W < ratio_activator_*k_activate_) {
            ratio_activator_ -= delta_catalyst_;
            ratio_deactivator_ += delta_catalyst_;
            atrpParticleP->updated = true;
            atrpParticleP->property_id = rc.property_id;
            atrpParticleP->p_state += rc.delta_state;
            atrpParticleP->p_type = prop->type();
            act++;
          }
        }
      }
    }
    for (std::map<longint, ATRPParticleP*>::iterator it = all_pids_state.begin(); it != all_pids_state.end(); it++) {
      if (it->second->updated) {
        selected_pids.push_back(it->second->p_id);
        selected_pids.push_back(it->second->p_state);
        selected_pids.push_back(it->second->property_id);
      }
      delete it->second;  // let's clean up
    }
    std::cout << "accept: " << accept << " reject: " << (total_N - accept) << " activate: " << act << " deactivate: " << dact << std::endl;
  } else {
    mpi::reduce(*(system.comm), my_N, total_N, std::plus<longint>(), 0);
    mpi::gather(*(system.comm), local_type_pids_state, global_type_pids, 0);
  }

  // Broadcast selected particles id to all CPUs.
  // TODO(jakub): instead of broadcast all selected pids, send only to CPUs that posses those particles.
  mpi::broadcast(*(system.comm), selected_pids, 0);

  // On every CPU, iterate over list of pids. If particle is on this CPU then process it.
  std::vector<Particle *> modified_particles;
  for (std::vector<longint>::iterator it = selected_pids.begin(); it != selected_pids.end();) {
    longint p_pid = *(it++);
    longint final_state = *(it++);
    longint property_id = *(it++);
    Particle *p = system.storage->lookupRealParticle(p_pid);  // check if particle is on this CPU.
    if (p) {
      shared_ptr<TopologyParticleProperties> pp = property_map_[property_id];
      p->state() = final_state;
      if (pp->updateParticleProperties(p)) {
        modified_particles.push_back(p);
      } else {
        throw std::runtime_error("Could not update properties of particle. Strange!");
      }
    }
  }
  // Synchronize all processess
  (*system.comm).barrier();

  timeUpdateParticles += wallTimer.getElapsedTime() - time0;

  // Update neighbour ghosts
  updateGhost(modified_particles);

  LOG4ESPP_INFO(theLogger, "ATRPActivator::updateParticles end");
}


/** Performs two-way parallel communication to update the ghost particles.
 * The parallel scheme is taken from
 * storage::DomainDecomposition::doGhostCommunication
 */
void ATRPActivator::updateGhost(const std::vector<Particle *> &modified_particles) {  // NOLINT
  real time0 = wallTimer.getElapsedTime();
  LOG4ESPP_DEBUG(theLogger, "ATRPActivator::updateGhost begin");

  int kCrCommTag = 0x69;

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
  real p_mass, p_q;

  out_buffer.write(data_length);

  for (std::vector<Particle *>::const_iterator it = modified_particles.begin(); it != modified_particles.end(); ++it) {
    p_id = (*it)->id();
    p_type = (*it)->type();
    p_state = (*it)->state();
    p_mass = (*it)->mass();
    p_q = (*it)->q();
    out_buffer.write(p_id);
    out_buffer.write(p_type);
    out_buffer.write(p_state);
    out_buffer.write(p_mass);
    out_buffer.write(p_q);
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
          in_buffer_0.read(p_mass);
          in_buffer_0.read(p_q);
        } else {
          in_buffer_1.read(p_id);
          in_buffer_1.read(p_type);
          in_buffer_1.read(p_state);
          in_buffer_1.read(p_mass);
          in_buffer_1.read(p_q);
        }

        // Update the ghost particle data on neighbour CPUs.
        particle = system.storage->lookupLocalParticle(p_id);

        if (particle && particle->ghost()) {
          LOG4ESPP_DEBUG(theLogger, "Update particle data");
          particle->setType(p_type);
          particle->setState(p_state);
          particle->setMass(p_mass);
          particle->setQ(p_q);
        }
      }
    }
    LOG4ESPP_DEBUG(theLogger, "Leaving unpack");
  }

  LOG4ESPP_DEBUG(theLogger, "ATRPActivator::updateGhost end");
  timeUpdateGhost += wallTimer.getElapsedTime() - time0;
}

void ATRPActivator::saveStatistics(std::string filename) {
  System &system = getSystemRef();
  if (system.comm->rank() == 0) {
    std::fstream fs(filename.c_str(), std::fstream::out);
    fs << "# step ratio_activator ratio_deactivator" << std::endl;
    for (std::vector<real>::iterator it = stats_k_activator.begin(); it != stats_k_activator.end();) {
      int step = (int) *(it++);
      real ratio_activator = *(it++);
      real ratio_deactivator = *(it++);
      fs << step << " " << std::scientific << ratio_activator << " " << std::scientific << ratio_deactivator
         << std::endl;
    }
    fs.close();
  }
}

std::vector<real> ATRPActivator::getStats() {
  std::vector<real> ret;
  ret.push_back(ratio_activator_);
  ret.push_back(ratio_deactivator_);
}

python::list ATRPActivator::getTimers() {
  python::list ret;

  ret.append(python::make_tuple("timeUpdateParticles", timeUpdateParticles));
  ret.append(python::make_tuple("timeUpdateGhost", timeUpdateGhost));
  ret.append(python::make_tuple("timeAll", timeUpdateGhost+timeUpdateParticles));

  return ret;
}

void ATRPActivator::registerPython() {
  using namespace espressopp::python;  // NOLINT

  class_<ATRPActivator, shared_ptr<ATRPActivator>, bases<Extension> >
      ("integrator_ATRPActivator", init<shared_ptr<System>, longint, longint, real, real, real, real, real>())
      .add_property("stats_filename", make_getter(&ATRPActivator::stats_filename_), make_setter(&ATRPActivator::stats_filename_))
      .add_property("select_from_all", make_getter(&ATRPActivator::select_from_all_), make_setter(&ATRPActivator::select_from_all_))
      .def("add_reactive_center", &ATRPActivator::addReactiveCenter)
      .def("update_particles", &ATRPActivator::updateParticles)
      .def("save_statistics", &ATRPActivator::saveStatistics)
      .def("connect", &ATRPActivator::connect)
      .def("disconnect", &ATRPActivator::disconnect);
}


}  // end namespace integrator
}  // end namespace espressopp
