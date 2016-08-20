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
#include "ChangeInRegion.hpp"

#include "types.hpp"
#include "System.hpp"
#include "storage/Storage.hpp"
#include "iterator/CellListIterator.hpp"

namespace espressopp {
using namespace iterator;  // NOLINT

namespace integrator {

LOG4ESPP_LOGGER(ChangeInRegion::theLogger, "ChangeInRegion");

ChangeInRegion::ChangeInRegion(shared_ptr<System> system, shared_ptr<ParticleRegion> particle_region)
    : Extension(system), particleRegion(particle_region) {
  LOG4ESPP_INFO(theLogger, "ChangeInRegion constructed");
}

void ChangeInRegion::disconnect() {
  sig_aftIntV.disconnect();
}

void ChangeInRegion::connect() {
  sig_aftIntV = integrator->aftIntV.connect(boost::bind(&ChangeInRegion::updateParticles, this));
}

void ChangeInRegion::updateParticles() {
  System& system = getSystemRef();
  for (ParticleRegion::iterator it=particleRegion->begin(); it != particleRegion->end(); it++) {
    Particle *p = *it;
    longint p_type = p->type();
    if (type_particleProperties.count(p_type) == 1) {
      type_particleProperties[p_type]->updateParticleProperties(p);
      LOG4ESPP_DEBUG(theLogger, "change property of particle " << p->id());
    }
    if (type_flags.count(p_type) == 1) {
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

/****************************************************
** REGISTRATION WITH PYTHON
****************************************************/
void ChangeInRegion::registerPython() {
  using namespace espressopp::python;  // NOLINT

  class_<ChangeInRegion, shared_ptr<ChangeInRegion>, bases<Extension> >
  ("integrator_ChangeInRegion", init<shared_ptr<System>, shared_ptr<ParticleRegion> >())
      .def("set_particle_properties", &ChangeInRegion::setParticleProperties)
      .def("set_flags", &ChangeInRegion::setFlags)
      .def("update_particles", &ChangeInRegion::updateParticles)
      .def("connect", &ChangeInRegion::connect)
      .def("disconnect", &ChangeInRegion::disconnect);
}

}  // end namespace integrator
}  // end namespace espressopp
