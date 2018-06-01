/*
  Copyright (C) 2018
      Zidan Zhang (zidan.zhang@kuleuven.be)
  
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
#include "GeneralReactionScheme.hpp"

#include "esutil/RNG.hpp"
#include "storage/Storage.hpp"


namespace espressopp {
namespace integrator {

LOG4ESPP_LOGGER(GeneralReactionScheme::theLogger, "GeneralReactionScheme");

GeneralReactionScheme::GeneralReactionScheme(shared_ptr<System> system) : Extension(system) {
  LOG4ESPP_INFO(theLogger, "GeneralReactionScheme constructed");


}

void GeneralReactionScheme::disconnect() {
  sig_aftIntV.disconnect();
}

void GeneralReactionScheme::connect() {
}

/**
 * The function get's the particles from CPU and do some operations
 */
void GeneralReactionScheme::updateParticles() {
  if (integrator->getStep() % (interval_) != 0)
    return;

  System &system = getSystemRef();

  // Collect on every CPU the particle properties.
  CellList cl = system.storage->getRealCells();
  longint my_N = 0;

  // TODO(jakub): this is expensive part
  for (espressopp::iterator::CellListIterator cit(cl); !cit.isDone(); ++cit) {
    Particle &p = *cit;
    my_N++;
  }

  // Send to master CPU
  longint my_rank = system.comm->rank();

  if (my_rank == 0) {

  } else {

  }

}


void GeneralReactionScheme::registerPython() {
  using namespace espressopp::python;  // NOLINT

  class_<GeneralReactionScheme, shared_ptr<GeneralReactionScheme>, bases<Extension> >
  ("integrator_GeneralReactionScheme", init<shared_ptr<System> >())
      .def("connect", &GeneralReactionScheme::connect)
      .def("disconnect", &GeneralReactionScheme::disconnect);
}

}  // end namespace integrator
}  // end namespace espressopp
