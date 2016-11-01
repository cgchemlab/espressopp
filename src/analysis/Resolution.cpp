/*
  Copyright (C) 2015
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
#include "Resolution.hpp"


using namespace espressopp;  //NOLINT
using namespace iterator;  // NOLINT

namespace espressopp {
namespace analysis {

real Resolution::compute_real() const {
  System &system = getSystemRef();
  CellList realCells = system.storage->getRealCells();

  real res = 0.0;
  longint NPart = 0;

  for (CellListIterator cit(realCells); !cit.isDone(); ++cit) {
    const Particle& p = *cit;
    res += p.lambda();
    NPart++;
  }

  longint systemN = 0;
  real sumRes = 0.0;

  mpi::all_reduce(*getSystem()->comm, res, sumRes, std::plus<real>());
  mpi::all_reduce(*getSystem()->comm, NPart, systemN, std::plus<longint>());

  return sumRes/systemN;
}

void Resolution::registerPython() {
  using namespace espressopp::python;  //NOLINT
  class_<Resolution, bases<Observable> >
    ("analysis_Resolution", init< shared_ptr<System> >())
    .add_property("value", &Resolution::compute_real);
}

real ResolutionFixedPairList::compute_real() const {
  real lambda = 0.0;
  longint Nbond = 0;

  for (FixedPairListLambda::IteratorParticleLambda it(fpl_->getParticlePairs()); it.isValid(); ++it) {
    lambda += it->lambda;
    Nbond++;
  }

  longint systemNbond = 0;
  real sumLambda = 0.0;

  mpi::all_reduce(*getSystem()->comm, lambda, sumLambda, std::plus<real>());
  mpi::all_reduce(*getSystem()->comm, Nbond, systemNbond, std::plus<longint>());

  return sumLambda/systemNbond;
}

void ResolutionFixedPairList::registerPython() {
  using namespace espressopp::python;  //NOLINT
  class_<ResolutionFixedPairList, bases<Observable> >
      ("analysis_ResolutionFixedPairList", init< shared_ptr<System>, shared_ptr<FixedPairListLambda> >())
      .add_property("value", &ResolutionFixedPairList::compute_real);
}

}  // end namespace analysis
}  // end namespace espressopp
