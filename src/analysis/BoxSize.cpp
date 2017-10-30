/*
  Copyright (C) 2017
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
#include "BoxSize.hpp"
#include "storage/DomainDecomposition.hpp"

using namespace espressopp;  // NOLINT

namespace espressopp {
namespace analysis {

std::vector<real> BoxSize::compute_real_vector() {

  Real3D box_sizeL = getSystem()->bc->getBoxL();

  return std::vector<real> {box_sizeL[0], box_sizeL[1], box_sizeL[2]};
}

void BoxSize::registerPython() {
  using namespace espressopp::python;  // NOLINT
  class_<BoxSize, bases<Observable> >
    ("analysis_BoxSize", init<shared_ptr<System> >());
}

}  // namespace analysis
}  // namespace espressopp
