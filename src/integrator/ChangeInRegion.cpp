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
using namespace iterator;

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
  for (ParticleRegion::iterator it=particleRegion->begin(); it != particleRegion->end(); it++ ) {
    Particle *p = *it;
    if (type_particleProperties.count(p->type()) == 1) {
      type_particleProperties[p->type()]->updateParticleProperties(p);
    }
  }
}

/****************************************************
** REGISTRATION WITH PYTHON
****************************************************/
void ChangeInRegion::registerPython() {

  using namespace espressopp::python;

  class_<ChangeInRegion, shared_ptr<ChangeInRegion>, bases<Extension> >
  ("integrator_ChangeInRegion", init<shared_ptr<System>, shared_ptr<ParticleRegion> >())
      .def("set_particle_properties", &ChangeInRegion::setParticleProperties)
      .def("update_particles", &ChangeInRegion::updateParticles)
      .def("connect", &ChangeInRegion::connect)
      .def("disconnect", &ChangeInRegion::disconnect);
}

}
}
