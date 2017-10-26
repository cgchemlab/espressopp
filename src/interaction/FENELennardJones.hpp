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
#ifndef _INTERACTION_FENELennardJones_HPP
#define _INTERACTION_FENELennardJones_HPP

#include <cmath>
#include "FixedPairListInteractionTemplate.hpp"
#include "Potential.hpp"

namespace espressopp {
namespace interaction {
/** This class provides methods to compute forces and energies of
    the FENELennardJones potential.

    \f[ V(r) = -\frac{1}{2} \Delta r_{max}^2 K \log \left[ 1 -
    \left(\frac{r-r_0}{\Delta r_{max}} \right)^2 \right]
    \f] + LJ(sigma, epsilon)
*/
class FENELennardJones : public PotentialTemplate<FENELennardJones> {
 private:
  real K;
  real r0;
  real rMax;
  real rMaxSqr;
  real sigma;
  real epsilon;
  real ff1, ff2;
  real ef1, ef2;
  real cutoffLj2_;

 public:
  static void registerPython();

  FENELennardJones() : K(0.0), r0(0.0), rMax(0.0), sigma(0.0), epsilon(0.0) {
    setShift(0.0);
    setCutoff(infinity);
    preset();
  }

  FENELennardJones(real _K, real _r0, real _rMax, real _sigma, real _epsilon, real _cutoff,
                   real _shift)
      : K(_K), r0(_r0), rMax(_rMax), sigma(_sigma), epsilon(_epsilon) {
    setShift(_shift);
    setCutoff(_cutoff);
    preset();
  }

  FENELennardJones(real _K, real _r0, real _rMax, real _sigma, real _epsilon, real _cutoff)
      : K(_K), r0(_r0), rMax(_rMax), sigma(_sigma), epsilon(_epsilon) {
    autoShift = false;
    setCutoff(_cutoff);
    setAutoShift();
    preset();
  }

  void preset() {
    rMaxSqr = rMax * rMax;
    real sig2 = sigma * sigma;
    real sig6 = sig2 * sig2 * sig2;
    ff1 = 48.0 * epsilon * sig6 * sig6;
    ff2 = 24.0 * epsilon * sig6;
    ef1 = 4.0 * epsilon * sig6 * sig6;
    ef2 = 4.0 * epsilon * sig6;
    cutoffLj2_ = std::pow(sigma, 1.0/3.0);
  }

  // Setter and getter
  void setK(real _K) {
    K = _K;
    updateAutoShift();
  }
  real getK() const { return K; }

  void setR0(real _r0) {
    r0 = _r0;
    updateAutoShift();
  }
  real getR0() const { return r0; }

  void setRMax(real _rMax) {
    rMax = _rMax;
    updateAutoShift();
    preset();
  }
  real getRMax() const { return rMax; }

  // Setter and getter
  void setEpsilon(real _epsilon) {
    epsilon = _epsilon;
    LOG4ESPP_INFO(theLogger, "epsilon=" << epsilon);
    updateAutoShift();
    preset();
  }

  real getEpsilon() const { return epsilon; }

  void setSigma(real _sigma) {
    sigma = _sigma;
    LOG4ESPP_INFO(theLogger, "sigma=" << sigma);
    updateAutoShift();
    preset();
  }
  real getSigma() const { return sigma; }

  real _computeEnergySqrRaw(real distSqr) const {
    real energy = -0.5 * rMaxSqr * K * log(1 - pow((sqrt(distSqr) - r0) / rMax, 2));
    if (distSqr <= cutoffLj2_) {
      real frac2 = sigma * sigma / distSqr;
      real frac6 = frac2 * frac2 * frac2;
      energy += 4.0 * epsilon * (frac6 * frac6 - frac6);
    }
    return energy;
  }

  bool _computeForceRaw(Real3D& force, const Real3D& dist, real distSqr) const {
    real ffactor;

    if (r0 > ROUND_ERROR_PREC) {
      real r = sqrt(distSqr);
      ffactor = -K * (r - r0) / (r * (1 - ((r - r0) * (r - r0) / rMaxSqr)));
    } else {
      ffactor = -K / (1.0 - distSqr / rMaxSqr);
    }
    force = dist * ffactor;

    // LJ force
    if (distSqr <= cutoffLj2_) {
      real frac2 = 1.0 / distSqr;
      real frac6 = frac2 * frac2 * frac2;
      real ffactorLJ = frac6 * (ff1 * frac6 - ff2) * frac2;
      force += dist * ffactorLJ;
    }

    return true;
  }

  static LOG4ESPP_DECL_LOGGER(theLogger);
};
}
}

#endif
