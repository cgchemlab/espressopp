/*
  Copyright (C) 2016-2017
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
#include "ChangeInRegion.hpp"

#include "esutil/RNG.hpp"
#include "storage/Storage.hpp"

#include <fstream>

namespace espressopp {
using namespace iterator;  // NOLINT

namespace integrator {

LOG4ESPP_LOGGER(ChangeInRegion::theLogger, "ChangeInRegion");

ChangeInRegion::ChangeInRegion(
    shared_ptr<System> system,
    shared_ptr<ParticleRegion> particle_region,
    real p)
      : Extension(system), particleRegion(particle_region) {
  LOG4ESPP_INFO(theLogger, "ChangeInRegion constructed");
  p_ = p;
  num_particles_ = 0;
  percentage_ = 0.0;
}

ChangeInRegion::ChangeInRegion(
    shared_ptr<System> system,
    shared_ptr<ParticleRegion> particle_region,
    long num_part)
    : Extension(system), particleRegion(particle_region) {
  LOG4ESPP_INFO(theLogger, "ChangeInRegion constructed");
  p_ = 0.0;
  num_particles_ = num_part;
  percentage_ = 0.0;
}

ChangeInRegion::ChangeInRegion(
    shared_ptr<System> system,
    shared_ptr<ParticleRegion> particle_region,
    long num_part, real percentage)
    : Extension(system), particleRegion(particle_region) {
  LOG4ESPP_INFO(theLogger, "ChangeInRegion constructed");
  p_ = 0.0;
  num_particles_ = 0;
  percentage_ = percentage;
}

void ChangeInRegion::disconnect() {
  sig_aftIntV.disconnect();
}

void ChangeInRegion::connect() {
  if (num_particles_ > 0)
    sig_aftIntV = integrator->aftIntV.connect(boost::bind(&ChangeInRegion::updateParticlesFlux, this));
  else if (p_ >  0.0)
    sig_aftIntV = integrator->aftIntV.connect(boost::bind(&ChangeInRegion::updateParticlesProb, this));
  else if (percentage_ > 0.0)
    sig_aftIntV = integrator->aftIntV.connect(boost::bind(&ChangeInRegion::updateParticlesPercentage, this));
}

void ChangeInRegion::updateParticlesFlux() {
  System& system = getSystemRef();
  shared_ptr< esutil::RNG > rng = system.rng;
  // collect particles on every cpu
  // send to master
  std::vector<longint> pids;
  for (ParticleRegion::iterator it=particleRegion->begin(); it != particleRegion->end(); it++) {
    Particle *p = *it;
    longint p_type = p->type();
    if (type_particleProperties.count(p_type) == 1) {
      pids.push_back(p->id());
    }
  }

  std::vector<std::vector<longint> > global_pids;
  std::vector<longint> selected_pids;
  if (system.comm->rank() == 0) {
    mpi::gather(*(system.comm), pids, global_pids, 0);

    std::vector<longint> all_pids;
    for (auto it = global_pids.begin(); it != global_pids.end(); ++it) {
      for (auto iit = it->begin(); iit != it->end(); ++iit) {
        all_pids.push_back(*iit);
      }
    }

    longint num_pids = all_pids.size();
    if (!stats_filename_.empty()) {
      std::fstream fs(stats_filename_.c_str(), std::fstream::out);
      fs << integrator->getStep() << " " << ((num_pids < num_particles_) ? num_pids : num_particles_) << std::endl;
      fs.close();
    }
    selected_pids.reserve(num_particles_);
    // select randomly num_particles_ pids;
    for (long i = 0; i < num_particles_ && num_pids > 0; i++) {
      longint idx = (*rng)(num_pids);
      selected_pids.push_back(all_pids[idx]);
      std::swap(all_pids[idx], all_pids[num_pids-1]);
      num_pids--;
    }
  } else {
    mpi::gather(*(system.comm), pids, global_pids, 0);
  }
  mpi::broadcast(*(system.comm), selected_pids, 0);

  for (auto it = selected_pids.begin(); it != selected_pids.end(); ++it) {
    Particle *p = system.storage->lookupRealParticle(*it);  // check if particle is on this CPU.
    if (p) {
      longint p_type = p->type();
      if (type_particleProperties.count(p_type) == 1) {
        type_particleProperties[p_type]->updateParticleProperties(p);
        LOG4ESPP_DEBUG(theLogger, "change property of particle " << p->id());
      }
      // action, depending on the type;
      if (type_flags.find(p_type) != type_flags.end()) {
        int flag = type_flags[p_type];
        LOG4ESPP_DEBUG(theLogger, "type_flags " << flag);
        if (flag & R_PARTICLE) {  // remove particle ;-)
          if (!p->ghost())
            system.storage->removeParticle(p->id());
        } else {
          if (flag & R_VELOCITY) {  // reset velocity
            p->setV(0.0);
            LOG4ESPP_DEBUG(theLogger, "reset velocity of particle " << p->id());
          }
          if (flag & R_FORCE) {  // reset force
            p->setF(0.0);
            LOG4ESPP_DEBUG(theLogger, "reset force of particle " << p->id());
          }
        }
      }
    }
  }
}

void ChangeInRegion::updateParticlesPercentage() {
  System& system = getSystemRef();
  shared_ptr< esutil::RNG > rng = system.rng;
  // collect particles on every cpu
  // send to master
  std::vector<longint> pids;
  for (ParticleRegion::iterator it=particleRegion->begin(); it != particleRegion->end(); it++) {
    Particle *p = *it;
    longint p_type = p->type();
    if (type_particleProperties.count(p_type) == 1) {
      pids.push_back(p->id());
    }
  }

  std::vector<std::vector<longint> > global_pids;
  std::vector<longint> selected_pids;
  if (system.comm->rank() == 0) {
    mpi::gather(*(system.comm), pids, global_pids, 0);

    std::vector<longint> all_pids;
    for (auto it = global_pids.begin(); it != global_pids.end(); ++it) {
      for (auto iit = it->begin(); iit != it->end(); ++iit) {
        all_pids.push_back(*iit);
      }
    }

    longint num_particles = all_pids.size();
    longint num_pids = round(percentage_*num_particles);
    selected_pids.reserve(num_pids);
    // select randomly num_particles_ pids;
    while (num_pids > 0) {
      longint idx = (*rng)(num_particles);
      selected_pids.push_back(all_pids[idx]);
      std::swap(all_pids[idx], all_pids[num_pids-1]);
      num_pids--;
    }
    if (!stats_filename_.empty()) {
      std::fstream fs(stats_filename_.c_str(), std::fstream::out);
      fs << integrator->getStep() << " " << percentage_*num_particles << std::endl;
      fs.close();
    }
  } else {
    mpi::gather(*(system.comm), pids, global_pids, 0);
  }
  mpi::broadcast(*(system.comm), selected_pids, 0);

  for (auto it = selected_pids.begin(); it != selected_pids.end(); ++it) {
    Particle *p = system.storage->lookupRealParticle(*it);  // check if particle is on this CPU.
    if (p) {
      longint p_type = p->type();
      if (type_particleProperties.count(p_type) == 1) {
        type_particleProperties[p_type]->updateParticleProperties(p);
        LOG4ESPP_DEBUG(theLogger, "change property of particle " << p->id());
      }
      // action, depending on the type;
      if (type_flags.find(p_type) != type_flags.end()) {
        int flag = type_flags[p_type];
        LOG4ESPP_DEBUG(theLogger, "type_flags " << flag);
        if (flag & R_PARTICLE) {  // remove particle ;-)
          if (!p->ghost())
            system.storage->removeParticle(p->id());
        } else {
          if (flag & R_VELOCITY) {  // reset velocity
            p->setV(0.0);
            LOG4ESPP_DEBUG(theLogger, "reset velocity of particle " << p->id());
          }
          if (flag & R_FORCE) {  // reset force
            p->setF(0.0);
            LOG4ESPP_DEBUG(theLogger, "reset force of particle " << p->id());
          }
        }
      }
    }
  }
}

void ChangeInRegion::updateParticlesProb() {
  System& system = getSystemRef();
  shared_ptr< esutil::RNG > rng = system.rng;
  for (ParticleRegion::iterator it=particleRegion->begin(); it != particleRegion->end(); it++) {
    if ((*rng)() < p_) {
      Particle *p = *it;
      longint p_type = p->type();
      if (type_particleProperties.count(p_type) == 1) {
        type_particleProperties[p_type]->updateParticleProperties(p);
        LOG4ESPP_DEBUG(theLogger, "change property of particle " << p->id());
      }
      // action, depending on the type;
      if (type_flags.find(p_type) != type_flags.end()) {
        int flag = type_flags[p_type];
        LOG4ESPP_DEBUG(theLogger, "type_flags " << flag);
        if (flag & R_PARTICLE) {  // remove particle ;-)
          if (!p->ghost())
            system.storage->removeParticle(p->id());
        } else {
          if (flag & R_VELOCITY) {  // reset velocity
            p->setV(0.0);
            LOG4ESPP_DEBUG(theLogger, "reset velocity of particle " << p->id());
          }
          if (flag & R_FORCE) {  // reset force
            p->setF(0.0);
            LOG4ESPP_DEBUG(theLogger, "reset force of particle " << p->id());
          }
        }
      }
    }
  }
}

/****************************************************
** REGISTRATION WITH PYTHON
****************************************************/
void ChangeInRegion::registerPython() {
  using namespace espressopp::python;  // NOLINT

  class_<ChangeInRegion, shared_ptr<ChangeInRegion>, bases<Extension> >
  ("integrator_ChangeInRegion", init<shared_ptr<System>, shared_ptr<ParticleRegion>, real>())
      .def(init<shared_ptr<System>, shared_ptr<ParticleRegion>, long>())
      .def(init<shared_ptr<System>, shared_ptr<ParticleRegion>, long, real>())
      .add_property("stats_filename",
                    make_getter(&ChangeInRegion::stats_filename_),
                    make_setter(&ChangeInRegion::stats_filename_))
      .def("set_particle_properties", &ChangeInRegion::setParticleProperties)
      .def("set_flags", &ChangeInRegion::setFlags)
      .def("update_particles_prob", &ChangeInRegion::updateParticlesProb)
      .def("update_particles_flux", &ChangeInRegion::updateParticlesFlux)
      .def("connect", &ChangeInRegion::connect)
      .def("disconnect", &ChangeInRegion::disconnect);
}

}  // end namespace integrator
}  // end namespace espressopp
