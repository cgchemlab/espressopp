/*
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

#include <math.h>
#include <bc/BC.hpp>

#include "python.hpp"
#include "TDforce.hpp"
#include "types.hpp"
#include "System.hpp"
#include "storage/Storage.hpp"
#include "iterator/CellListIterator.hpp"
#include "interaction/InterpolationLinear.hpp"
#include "interaction/InterpolationAkima.hpp"
#include "interaction/InterpolationCubic.hpp"


namespace espresso {
using namespace iterator;  // NOLINT
namespace integrator {

LOG4ESPP_LOGGER(TDforce::theLogger, "TDforce");

TDforce::TDforce(shared_ptr<System> system) : Extension(system) {
  type = Extension::FreeEnergyCompensation;

  center = (0.0, 0.0, 0.0);
  LOG4ESPP_INFO(theLogger, "TDforce constructed");
}

TDforce::~TDforce() {}

void TDforce::connect() {
  _applyForce = integrator->aftCalcF.connect(
  boost::bind(&TDforce::applyForce, this));
}

void TDforce::disconnect() {
  _applyForce.disconnect();
}


void TDforce::addForce(int itype, const char* _filename, int type, int force_type) {
  boost::mpi::communicator world;
  filename = _filename;
  Table table;

  if (force_type != 0 && force_type != 1)
    throw std::runtime_error("Wrong TD force type. Set 0 for x-direction and 1 for radius");
  force_type_ = force_type;

  switch (itype) {
    case 1:  // InterpolationLinear
      table = make_shared <interaction::InterpolationLinear> ();
      table->read(world, _filename);
      break;
    case 2:  // InterpolationAkima
      table = make_shared <interaction::InterpolationAkima> ();
      table->read(world, _filename);
      break;
    case 3:  // InterpolationCubic
      table = make_shared <interaction::InterpolationCubic> ();
      table->read(world, _filename);
      break;
    default:
      throw std::runtime_error("Wrong type of interpolation.");
  }
  forces.insert(std::make_pair(type, table));
}


void TDforce::applyForce() {
  LOG4ESPP_DEBUG(theLogger, "apply TD force");

  System& system = getSystemRef();

  // iterate over CG particles
  CellList cells = system.storage->getRealCells();
  for (CellListIterator cit(cells); !cit.isDone(); ++cit) {
    Table table = forces.find(cit->getType())->second;
    if (table) {
      // calculate distance from reference point
      Real3D dl;
      system.bc->getMinimumImageVector(dl, cit->getPos(), center);

      if (force_type_ == 1) {
        real dist = dl.abs();

        // read fforce from table
        real fforce = table->getForce(dist);
        fforce /= dist;

        // substract td force
        cit->force() += (dl * fforce);

      } else if (force_type_ == 0) {
        // use this if you need 1-dir force only!
        real force = table->getForce(fabs(dl[0])) * (dl[0] / fabs(dl[0]));
        cit->force()[0] += force;
      }
    }
  }
}

void TDforce::setCenter(real x, real y, real z) {
  center = Real3D(x, y, z);
}

/****************************************************
** REGISTRATION WITH PYTHON
****************************************************/
void TDforce::registerPython() {
  using namespace espresso::python;  //NOLINT

  void (TDforce::*pySetCenter)(real x, real y, real z) = &TDforce::setCenter;

  void (TDforce::*pyAddForce)(int itype, const char* filename, int type, int force_type)
      = &TDforce::addForce;

  class_<TDforce, shared_ptr<TDforce>, bases<Extension> >
      ("integrator_TDforce", init< shared_ptr<System> >())
        .add_property("filename", &TDforce::getFilename)
        .def("connect", &TDforce::connect)
        .def("disconnect", &TDforce::disconnect)
        .def("setCenter", pySetCenter)
        .def("addForce", pyAddForce);
}
}  // namespace integrator
}  // namespace espresso
