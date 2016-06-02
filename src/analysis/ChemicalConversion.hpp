/*
  Copyright (C) 2016
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
#ifndef _ANALYSIS_CHEMICALCONVERSION_HPP
#define _ANALYSIS_CHEMICALCONVERSION_HPP

#include "types.hpp"
#include "Observable.hpp"
#include "ParticleGroup.hpp"
#include "storage/DomainDecomposition.hpp"
#include "iterator/CellListIterator.hpp"
#include "boost/signals2.hpp"
#include "boost/serialization/map.hpp"

namespace espressopp {
namespace analysis {

class ChemicalConversion : public Observable {
 public:
  ChemicalConversion(shared_ptr<System> system, longint _p_type, longint total) :
      Observable(system), total_value(total), p_type(_p_type) {
    result_type = real_scalar;
  }

  ~ChemicalConversion() {}
  real compute_real() const;

  boost::signals2::signal1<void, real> onValue;

  static void registerPython();
 private:
  real total_value;
  longint p_type;
};


class ChemicalConversionTypeSequence : public Observable {
 public:
  ChemicalConversionTypeSequence(shared_ptr<System> system, shared_ptr<ParticleGroup> pg_, longint total) :
      Observable(system), total_value(total), particle_group_(pg_) {
    result_type = real_scalar;
  }

  ~ChemicalConversionTypeSequence() {}
  real compute_real() const;

  void setSequence(std::vector<longint> in_seq) {
    type_seq_ = in_seq;
  }

  boost::signals2::signal1<void, real> onValue;

  static void registerPython();
 private:
  real total_value;
  std::vector<longint> type_seq_;
  shared_ptr<ParticleGroup> particle_group_;
};

}  // end namespace analysis
}  // end namespace espressopp

#endif
