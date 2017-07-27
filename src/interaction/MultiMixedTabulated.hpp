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
#ifndef _INTERACTION_MultiMixedTabulated_HPP
#define _INTERACTION_MultiMixedTabulated_HPP

#include <string>

#include "Potential.hpp"
#include "Interpolation.hpp"
#include "analysis/ChemicalConversion.hpp"

namespace espressopp {

namespace interaction {

class MultiMixedTabulated: public PotentialTemplate<MultiMixedTabulated> {
 public:
  static LOG4ESPP_DECL_LOGGER(theLogger);

  MultiMixedTabulated() {
    setShift(0.0);
    setCutoff(infinity);
    interpolationType = 0;
    mix_value_ = 0.0;
  }

  MultiMixedTabulated(int itype, real cutoff) {
    setShift(0.0);
    setCutoff(cutoff);
    interpolationType = itype;
    mix_value_ = 0.0;
  }

  void registerTableOnRange(
      std::string tab1,
      std::string tab2,
      shared_ptr<analysis::ChemicalConversion> chm,
      real min_value,
      real max_value);

  real _computeEnergySqrRaw(real distSqr) const {
    if (interpolationType != 0) {
      real dist = sqrt(distSqr);

      real e1 = current_table1->getEnergy(dist);
      real e2 = current_table2->getEnergy(dist);

      real val = mix_value_ * e1 + (1.0 - mix_value_) * e2;

      return val;
    } else {
      LOG4ESPP_ERROR(theLogger, "interpolationType = " << interpolationType);
      return 0.0;
    }
  }

  bool _computeForceRaw(Real3D &force, const Real3D &dist, real distSqr) const {
    real ffactor1, ffactor2, ffactor;
    if (interpolationType != 0) {
      real dist = sqrt(distSqr);
      ffactor1 = current_table1->getForce(dist);
      ffactor2 = current_table2->getForce(dist);

      ffactor = mix_value_*ffactor1 + (1.0 - mix_value_) * ffactor2;

      force = dist * ffactor;
      return true;
    } else {
      return false;
    }
  }

  static void registerPython();

 private:
  shared_ptr<Interpolation> current_table1;
  shared_ptr<Interpolation> current_table2;
  void onValue(
      real min_value,
      real max_value,
      shared_ptr<Interpolation> table1,
      shared_ptr<Interpolation> table2,
      real value);
  int interpolationType;
  real mix_value_;
};

}  // namespace interaction
}  // namespace espressopp
#endif
