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
#ifndef _ANALYSIS_ATRPACTIVATORSTATS_HPP
#define _ANALYSIS_ATRPACTIVATORSTATS_HPP

#include <integrator/ATRPActivator.hpp>
#include "types.hpp"
#include "Observable.hpp"
#include "ParticleGroup.hpp"
#include "storage/DomainDecomposition.hpp"
#include "iterator/CellListIterator.hpp"
#include "boost/signals2.hpp"
#include "boost/serialization/map.hpp"

namespace espressopp {
namespace analysis {

class ATRPActivatorStats : public Observable {
 public:
  ATRPActivatorStats(shared_ptr<System> system, shared_ptr<integrator::ATRPActivator> atrp_activator) :
      Observable(system), atrp_activator_(atrp_activator) {
    result_type = real_vector;
    result_vector_size = 2;
  }

  ~ATRPActivatorStats() {}
  real compute_real() const;
  std::vector<real> compute_real_vector();

  boost::signals2::signal<void (real)> onValue;

  static void registerPython();
 private:
  shared_ptr<integrator::ATRPActivator> atrp_activator_;
};

}  // end namespace analysis
}  // end namespace espressopp

#endif
