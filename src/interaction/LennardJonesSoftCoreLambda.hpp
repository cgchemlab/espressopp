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

#ifndef _INTERACTION_LENNARDJONESSOFTCORELAMBDA_HPP
#define _INTERACTION_LENNARDJONESSOFTCORELAMBDA_HPP

#include "FixedPairListInteractionTemplate.hpp"
#include "FixedPairListTypesInteractionTemplate.hpp"
#include "Potential.hpp"

namespace espressopp {
namespace interaction {
/** This class provides methods to compute forces and energies of
the Lennard Jones potential.

\f[ V(r) = 4 \varepsilon \left[ \left( \frac{\sigma}{r} \right)^{12} -
\left( \frac{\sigma}{r} \right)^{6} \right]
\f]

This class implement potential defined in:
Peters, J. H. et al. Phys. Rev. E - Stat. Nonlinear, Soft Matter Phys. 94, 1-19 (2016).

*/
class LennardJonesSoftCoreLambda: public PotentialTemplate<LennardJonesSoftCoreLambda> {
 private:
  real epsilon;
  real sigma;
  real sigma6_;
  real alpha_;
  real weight_power_;
  real ff1, ff2;
  real ef1, ef2;

 public:
  static void registerPython();

  LennardJonesSoftCoreLambda() : epsilon(0.0), sigma(0.0), alpha_(0.0) {
    setShift(0.0);
    setCutoff(infinity);
    preset();
  }

  LennardJonesSoftCoreLambda(real _epsilon, real _sigma, real _cutoff, real _shift, real alpha, real weight_power)
      : epsilon(_epsilon), sigma(_sigma), alpha_(alpha), weight_power_(weight_power) {
    setShift(_shift);
    setCutoff(_cutoff);
    preset();
  }

  LennardJonesSoftCoreLambda(real _epsilon, real _sigma, real _cutoff, real alpha, real weight_power)
      : epsilon(_epsilon), sigma(_sigma), alpha_(alpha), weight_power_(weight_power) {
    autoShift = false;
    setCutoff(_cutoff);
    preset();
    setAutoShift();
  }

  virtual ~LennardJonesSoftCoreLambda() {};

  void preset() {
    real sig2 = sigma * sigma;
    real sigma6_ = sig2 * sig2 * sig2;
    ff1 = 48.0 * epsilon * sigma6_ * sigma6_;
    ff2 = 24.0 * epsilon * sigma6_;
    ef1 =  4.0 * epsilon * sigma6_ * sigma6_;
    ef2 =  4.0 * epsilon * sigma6_;
  }

  // Setter and getter
  void setEpsilon(real _epsilon) {
    epsilon = _epsilon;
    LOG4ESPP_INFO(theLogger, "epsilon=" << epsilon);
    preset();
    updateAutoShift();
  }

  real getEpsilon() const { return epsilon; }

  void setSigma(real _sigma) {
    sigma = _sigma;
    sigma6_ = pow(sigma, 6);
    LOG4ESPP_INFO(theLogger, "sigma=" << sigma);
    preset();
    updateAutoShift();
  }
  real getSigma() const { return sigma; }

  real getAlpha() const { return alpha_; }
  void setAlpha(real s_) { alpha_ = s_;}
  real getWeightPower() const { return weight_power_; }
  void setWeightPower(real s_) { weight_power_ = s_; }


  real _computeEnergy(const Particle &p1, const Particle &p2) const {
    real p1_lambda = p1.lambda();
    real p2_lambda = p2.lambda();
    real lambda_sqr = (1.0 - p1_lambda * p2_lambda);

    Real3D dist = p1.position() - p2.position();
    real distSqr = dist.sqr();
    real dist6 = distSqr*distSqr*distSqr;

    real distEffSqr = pow(alpha_*sigma6_*pow(lambda_sqr, weight_power_) + dist6, 1.0/3.0);

    if (distEffSqr > cutoffSqr)
      return 0.0;

    real energy = _computeEnergySqrRaw(distEffSqr);
    return energy - shift;
  }

  real _computeEnergySqrRaw(real distSqr) const {
    real frac2 = sigma*sigma / distSqr;
    real frac6 = frac2 * frac2 * frac2;
    real energy = 4.0 * epsilon * (frac6 * frac6 - frac6);
    return energy;
  }

  bool _computeForce(Real3D &force, const Particle &p1, const Particle &p2) const {
    real p1_lambda = p1.lambda();
    real p2_lambda = p2.lambda();
    real lambda_sqr = (1.0 - p1_lambda * p2_lambda);

    Real3D dist = p1.position() - p2.position();
    Real3D e_dist = dist / dist.abs();

    real distSqr = dist.sqr();
    real dist6 = distSqr*distSqr*distSqr;

    real distEffAbs = pow(alpha_*sigma6_*pow(lambda_sqr, weight_power_) + dist6, 1.0/6.0);
    Real3D distEff = distEffAbs * e_dist;

    return _computeForceRaw(force, distEff, distEffAbs*distEffAbs);
  }


  bool _computeForceRaw(Real3D &force,
                        const Real3D &dist,
                        real distSqr) const {
    real frac2 = 1.0 / distSqr;
    real frac6 = frac2 * frac2 * frac2;
    real ffactor = frac6 * (ff1 * frac6 - ff2) * frac2;
    force = dist * ffactor;
    return true;;
  }

  static LOG4ESPP_DECL_LOGGER(theLogger);
};

// provide pickle support
struct LennardJonesSoftCoreLambda_pickle: boost::python::pickle_suite {
  static boost::python::tuple getinitargs(LennardJonesSoftCoreLambda const &pot) {
    real eps = pot.getEpsilon();
    real sig = pot.getSigma();
    real rc = pot.getCutoff();
    real sh = pot.getShift();
    real alpha = pot.getAlpha();
    real weight_power = pot.getWeightPower();
    return boost::python::make_tuple(eps, sig, rc, sh, alpha, weight_power);
  }
};
}
}

#endif
