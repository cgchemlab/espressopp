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
#include "NFixedListEntries.hpp"
#include "FixedPairList.hpp"

using namespace espressopp;  //NOLINT

namespace espressopp {
namespace analysis {

real NFixedPairListEntries::compute_real() const {
  return fpl_->totalSize();
}

void NFixedPairListEntries::registerPython() {
  using namespace espressopp::python;  //NOLINT
  class_<NFixedPairListEntries, bases<Observable> >
    ("analysis_NFixedPairListEntries",
        init< shared_ptr<System>, shared_ptr<FixedPairList> >())
    .add_property("value", &NFixedPairListEntries::compute_real);
}
}  // end namespace analysis
}  // end namespace espressopp
