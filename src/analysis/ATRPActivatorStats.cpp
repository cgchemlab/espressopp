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
#include "ATRPActivatorStats.hpp"

using namespace espressopp;  //NOLINT

namespace espressopp {
namespace analysis {

real ATRPActivatorStats::compute_real() const {
  throw std::runtime_error("This observable returns vector, use compute_real_vector.");
}

void ATRPActivatorStats::registerPython() {
  using namespace espressopp::python;  //NOLINT
  class_<ATRPActivatorStats, bases<Observable>, boost::noncopyable>
    ("analysis_ATRPActivatorStats",
        init< shared_ptr<System>, shared_ptr<integrator::ATRPActivator> >());
}
std::vector<real> ATRPActivatorStats::compute_real_vector() {
  result_real_vector = atrp_activator_->getStats();
  return result_real_vector;
}

}  // end namespace analysis
}  // end namespace espressopp
