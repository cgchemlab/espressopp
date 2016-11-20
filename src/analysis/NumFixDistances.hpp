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

// ESPP_CLASS
#ifndef _ANALYSIS_NFIXDISTANCES_HPP
#define _ANALYSIS_NFIXDISTANCES_HPP

#include "types.hpp"
#include "Observable.hpp"
#include "integrator/FixDistances.hpp"

namespace espressopp {
namespace analysis {

class NumFixDistances : public Observable {
 public:
  NumFixDistances(shared_ptr<System> system, shared_ptr<integrator::FixDistances> fd) :
      Observable(system), fd_(fd) {
    result_type = real_scalar;
  }

  ~NumFixDistances() {}
  real compute_real() const {
    return fd_->totalSize();
  }

  static void registerPython();
 private:
  shared_ptr<integrator::FixDistances> fd_;
};

}  // end namespace analysis
}  // end namespace espressopp

#endif
