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
#ifndef _INTERACTION_TABULATEDCAPPED_HPP
#define _INTERACTION_TABULATEDCAPPED_HPP

#include "Potential.hpp"
#include "Interpolation.hpp"
#include "../include/esconfig.hpp"

namespace espressopp {

namespace interaction {

/** This class provides methods to compute forces and energies of
a tabulated potential.

    The potential and forces must be provided in a file.

    Be careful: default and copy constructor of this class are used.
*/

class TabulatedCapped: public PotentialTemplate<TabulatedCapped> {

 private:
  std::string filename;
  shared_ptr <Interpolation> table;
  int interpolationType;
  real caprad_;
  real capradSqr_;

  real caprad_energy_;
  real caprad_force_;

 public:
  static void registerPython();

  TabulatedCapped() {
    setShift(0.0);
    setCutoff(infinity);
    interpolationType = 0;
    caprad_ = 0.0;
    capradSqr_ = 0.0;
    initialized = false;
  }

  // used for fixedpairlist (2-body bonded interaction)
  TabulatedCapped(int itype, const char *filename) {
    setInterpolationType(itype);
    setFilename(itype, filename);
    setShift(0.0);
    setCutoff(infinity);
    caprad_ = 0.0;
    capradSqr_ = 0.0;
    initialized = true;

    preset();
  }

  TabulatedCapped(int itype, const char *filename, real cutoff) {
    setInterpolationType(itype);
    setFilename(itype, filename);
    setShift(0.0);
    setCutoff(cutoff);
    caprad_ = 0.0;
    capradSqr_ = 0.0;
    initialized = true;

    preset();
  }

  TabulatedCapped(int itype, const char *filename, real cutoff, real caprad) {
    setInterpolationType(itype);
    setFilename(itype, filename);
    setShift(0.0);
    setCutoff(cutoff);
    caprad_ = caprad;
    capradSqr_ = caprad * caprad;
    initialized = true;

    preset();
  }

  void preset() {
    caprad_energy_ = table->getEnergy(caprad_);
    caprad_force_ = table->getForce(caprad_);
  }

  real getCaprad() const;

  /** Setter for the interpolation type */
  void setInterpolationType(int itype) { interpolationType = itype; }

  /** Getter for the interpolation type */
  int getInterpolationType() const { return interpolationType; }

  /** Setter for the filename will read in the table. */
  void setFilename(int itype, const char *_filename);

  /** Getter for the filename. */
  const char *getFilename() const { return filename.c_str(); }

  real _computeEnergySqrRaw(real distSqr) const {
    // make an interpolation
    if (interpolationType != 0) {
      if (distSqr <= capradSqr_)
        return caprad_energy_;
      return table->getEnergy(sqrt(distSqr));
    } else {
      return 0;
    }
  }

  bool _computeForceRaw(Real3D &force, const Real3D &dist, real distSqr) const {
    real ffactor = 0.0;
    if (interpolationType != 0) {
      if (distSqr < capradSqr_) {
        force = caprad_force_ * (dist/dist.abs());
        return true;
      } else {
        real distrt = sqrt(distSqr);
        ffactor = table->getForce(distrt);
        ffactor /= distrt;
      }
    } else {
      return false;
    }
    force = dist * ffactor;
    return true;
  }

};//class

// provide pickle support
struct TabulatedCapped_pickle: boost::python::pickle_suite {
  static
  boost::python::tuple
  getinitargs(TabulatedCapped const &pot) {
    int itp = pot.getInterpolationType();
    std::string fn = pot.getFilename();
    real rc = pot.getCutoff();
    real caprad = pot.getCaprad();
    return boost::python::make_tuple(itp, fn, rc, caprad);
  }
};

}
}

#endif
