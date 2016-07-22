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
#include "MaxAvgMinForce.hpp"
#include "storage/DomainDecomposition.hpp"
#include "iterator/CellListIterator.hpp"

using namespace espressopp;
using namespace iterator;

namespace espressopp {
namespace analysis {

real MaxForce::compute() const {
  System &system = getSystemRef();
  CellList realCells = system.storage->getRealCells();

  real maxForce = 0.0;
  for (CellListIterator cit(realCells); !cit.isDone(); ++cit) {
    real forceAbs = cit->force().abs();
    if (forceAbs > maxForce && !(isnan(forceAbs) || isinf(forceAbs))) {
      maxForce = forceAbs;
    }
  }

  // It was reduce.
  real systemMaxForce = 0.0;
  boost::mpi::all_reduce(*getSystem()->comm, maxForce, systemMaxForce, mpi::maximum<real>());

  return systemMaxForce;
}

void MaxForce::registerPython() {
  using namespace espressopp::python;
  class_<MaxForce, bases<Observable> >
    ("analysis_MaxForce", init<shared_ptr<System> >());
}

}
}
