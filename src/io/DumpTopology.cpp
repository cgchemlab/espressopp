/*
  Copyright (c) 2015
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
#include <fstream>
#include <iomanip>
#include "DumpPairs.hpp"
#include "storage/Storage.hpp"

#include "bc/BC.hpp"

#include "analysis/ConfigurationExt.hpp"
#include "analysis/ConfigurationsExt.hpp"

using namespace espressopp;
using namespace espressopp::analysis;
using namespace std;

namespace espressopp {
namespace io {

void DumpPairs::dump() {
  shared_ptr <System> system = getSystem();

}

// Python wrapping
void DumpPairs::registerPython() {

  using namespace espressopp::python;

  class_ < DumpPairs, bases < ParticleAccess >, boost::noncopyable >
      ("io_DumpPairs", init < shared_ptr < System >,
          shared_ptr < integrator::MDIntegrator >,
          std::string,
          bool,
          real,
          std::string,
          bool > ())
          .add_property("filename", &DumpPairs::getFilename,
                        &DumpPairs::setFilename)
          .add_property("unfolded", &DumpPairs::getUnfolded,
                        &DumpPairs::setUnfolded)
          .add_property("append", &DumpPairs::getAppend,
                        &DumpPairs::setAppend)
          .add_property("length_factor", &DumpPairs::getLengthFactor,
                        &DumpPairs::setLengthFactor)
          .add_property("length_unit", &DumpPairs::getLengthUnit,
                        &DumpPairs::setLengthUnit)
          .def("dump", &DumpPairs::dump);
}
}
}
