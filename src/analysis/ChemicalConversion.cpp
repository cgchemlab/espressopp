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
#include "ChemicalConversion.hpp"

using namespace espressopp;  //NOLINT

namespace espressopp {
namespace analysis {

real ChemicalConversion::compute_real() const {
  System& system = getSystemRef();
  CellList realCells = system.storage->getRealCells();

  longint local_count = 0;
  for (iterator::CellListIterator cit(realCells); !cit.isDone(); ++cit) {
    if (cit->type() == p_type)
      local_count++;
  }

  longint global_count = 0;
  boost::mpi::all_reduce(*getSystem()->comm, local_count, global_count, std::plus<longint>());

  real result = global_count / total_value;
  // Send value via signal.
  onValue(result);

  return result;
}

void ChemicalConversion::registerPython() {
  using namespace espressopp::python;  //NOLINT
  class_<ChemicalConversion, bases<Observable>, boost::noncopyable>
    ("analysis_ChemicalConversion",
        init< shared_ptr<System>, longint, longint >())
    .add_property("value", &ChemicalConversion::compute_real);
}
}  // end namespace analysis
}  // end namespace espressopp
