/*
  Copyright (C) 20126
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
#ifndef _INTERACTION_MultiTabulated_HPP
#define _INTERACTION_MultiTabulated_HPP

#include "Potential.hpp"
#include "Interpolation.hpp"
#include "analysis/ChemicalConversion.hpp"

namespace espressopp {

namespace interaction {

/** This class provides methods to compute forces and energies of a Tabulated potential.
    The potential and forces must be provided in a file.

    Be careful: default and copy constructor of this class are used.
*/

class MultiTabulated: public PotentialTemplate<MultiTabulated> {
 public:
  static LOG4ESPP_DECL_LOGGER(theLogger);

  MultiTabulated() {
    setShift(0.0);
    setCutoff(infinity);
  }

  MultiTabulated(real cutoff) {
    setShift(0.0);
    setCutoff(cutoff);
  }

  void registerTableOnRange(
      const char *filename,
      int itype,
      shared_ptr<analysis::ChemicalConversion> chm,
      real min_value,
      real max_value,
      bool default_table
  );

  real _computeEnergySqrRaw(real distSqr) const {
    if (current_table)
      return current_table->getEnergy(sqrt(distSqr));
    else
      return 0;
  }

  bool _computeForceRaw(Real3D &force, const Real3D &dist, real distSqr) const {
    real ffactor;
    if (current_table) {
      real distrt = sqrt(distSqr);
      ffactor = current_table->getForce(distrt);
      ffactor /= distrt;
    } else {
      return false;
    }
    force = dist * ffactor;
    return true;
  }

  static void registerPython();

 private:
  shared_ptr<Interpolation> current_table;
  void onValue(real min_value, real max_value, shared_ptr<Interpolation> table, real value);
};

}
}
#endif
