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
#ifndef _ANALYSIS_NFIXEDLISTENTRIES_HPP
#define _ANALYSIS_NFIXEDLISTENTRIES_HPP

#include "types.hpp"
#include "Observable.hpp"
#include "FixedPairList.hpp"
#include "FixedTripleList.hpp"
#include "FixedQuadrupleList.hpp"

namespace espressopp {
namespace analysis {

class NFixedListEntries : public Observable {
 public:
  NFixedListEntries(shared_ptr<System> system, shared_ptr<FixedPairList> fpl) :
      Observable(system), fpl_(fpl) {
    result_type = real_scalar;
  }
  NFixedListEntries(shared_ptr<System> system, shared_ptr<FixedTripleList> ftl) :
      Observable(system), ftl_(ftl) {
    result_type = real_scalar;
  }
  NFixedListEntries(shared_ptr<System> system, shared_ptr<FixedQuadrupleList> fql) :
      Observable(system), fql_(fql) {
    result_type = real_scalar;
  }

  ~NFixedListEntries() {}
  real compute_real() const;

  static void registerPython();
 private:
  shared_ptr<FixedPairList> fpl_;
  shared_ptr<FixedTripleList> ftl_;
  shared_ptr<FixedQuadrupleList> fql_;
};
}  // end namespace analysis
}  // end namespace espressopp

#endif
