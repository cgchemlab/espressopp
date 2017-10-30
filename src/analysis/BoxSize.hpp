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

// ESPP_CLASS
#ifndef _ANALYSIS_BOXSIZE_HPP
#define _ANALYSIS_BOXSIZE_HPP

#include "types.hpp"
#include "Observable.hpp"
#include "bc/BC.hpp"

namespace espressopp {
namespace analysis {

/** Class to get the maximum force acting on particles in the system. */
class BoxSize : public Observable {
 public:
  explicit BoxSize(shared_ptr<System> system) : Observable(system) {
    result_type = real_vector;
    result_vector_size = 3;
  }
  virtual ~BoxSize() {}
  real compute_real() const {
    throw std::runtime_error("This observable returns vector, use compute_real_vector.");
  }
  std::vector<real> compute_real_vector();

  static void registerPython();
};

}  // namespace analysis
}  // namespace espressopp
#endif
