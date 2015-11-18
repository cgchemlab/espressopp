/*
  Copyright (C) 2015
      Jakub Krajniak (jkrajniak at gmail.com)
  Copyright (C) 2012,2013
      Max Planck Institute for Polymer Research
  Copyright (C) 2008,2009,2010,2011
      Max-Planck-Institute for Polymer Research & Fraunhofer SCAI
  
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
#ifndef _INTERACTION_LENNARDJONESLAMBDA_HPP
#define _INTERACTION_LENNARDJONESLAMBDA_HPP

#include "FixedPairListInteractionTemplate.hpp"
#include "FixedPairListTypesInteractionTemplate.hpp"
#include "Potential.hpp"

namespace espressopp {
  namespace interaction {
    /** This class provides methods to compute forces and energies of
	the Lennard Jones potential.

	\f[ V(r) = 4 \varepsilon \left[ \left( \frac{\lambda \sigma}{r} \right)^{12} -
	\left( \frac{\lambda \sigma}{r} \right)^{6} \right]
	\f]

    */
    class LennardJonesLambda : public PotentialTemplate< LennardJonesLambda > {
    private:
      real epsilon;
      real sigma;
      real max_force;
      bool has_max_force;

    public:
      static void registerPython();

      LennardJonesLambda() : epsilon(0.0), sigma(0.0) {
        setShift(0.0);
        setCutoff(infinity);
      }

      LennardJonesLambda(real _epsilon, real _sigma, real _cutoff, real _shift)
          : epsilon(_epsilon), sigma(_sigma) {
        setShift(_shift);
        setCutoff(_cutoff);
        has_max_force = false;
      }

      LennardJonesLambda(real _epsilon, real _sigma, real _cutoff)
          : epsilon(_epsilon), sigma(_sigma) {
        autoShift = false;
        setCutoff(_cutoff);
        setAutoShift();
        has_max_force = false;
      }

      virtual ~LennardJonesLambda() {};

      // Setter and getter
      void setEpsilon(real _epsilon) {
        epsilon = _epsilon;
        LOG4ESPP_INFO(theLogger, "epsilon=" << epsilon);
        updateAutoShift();
      }
      
      real getEpsilon() const { return epsilon; }

      void setSigma(real _sigma) { 
        sigma = _sigma; 
        LOG4ESPP_INFO(theLogger, "sigma=" << sigma);
        updateAutoShift();
      }
      real getSigma() const { return sigma; }

      real getMaxForce() const { return max_force; }
      void setMaxForce(real _maxForce) {
        max_force = _maxForce;
        has_max_force = _maxForce != -1;
      }

      real _computeEnergy(const Particle& p1, const Particle& p2) const {
        real p1_lambda = p1.lambda();
        real p2_lambda = p2.lambda();
        real lambda_sqr = p1_lambda * p2_lambda;
        if (lambda_sqr == 0.0)
          return 0.0;

        Real3D dist = p1.position() - p2.position();
        real distSqr = dist.sqr();
        if (distSqr > cutoffSqr)
          return 0.0;
        real frac2 = sigma*sigma*lambda_sqr / dist.sqr();
        real frac6 = frac2 * frac2 * frac2;
        real energy = 4.0 * epsilon * (frac6 * frac6 - frac6);
        return energy - shift;
      }

      real _computeEnergySqrRaw(real distSqr) const {
        real frac2 = sigma*sigma / distSqr;
        real frac6 = frac2 * frac2 * frac2;
        real energy = 4.0 * epsilon * (frac6 * frac6 - frac6);
        return energy;
      }

      bool _computeForce(Real3D& force, const Particle &p1, const Particle &p2) const {
        real p1_lambda = p1.lambda();
        real p2_lambda = p2.lambda();
        real lambda_sqr = p1_lambda * p2_lambda;
        if (lambda_sqr == 0.0)
          return false;

        Real3D dist = p1.position() - p2.position();
        real sig2 = sigma * sigma * lambda_sqr;
        real sig6 = sig2 * sig2 * sig2;
        real ff1_ = 48.0 * epsilon * sig6 * sig6;
        real ff2_ = 24.0 * epsilon * sig6;

        real frac2 = 1.0 / dist.sqr();
        real frac6 = frac2 * frac2 * frac2;
        real ffactor = frac6 * (ff1_ * frac6 - ff2_) * frac2;
        force = dist * ffactor;

        if (lambda_sqr != 1.0 && has_max_force) {
          if (force.isNaNInf()) {  // Force is inf.
            force = dist * max_force;
          } else {
            real abs_force = force.abs();
            if (abs_force > max_force) {
              force = dist * max_force;
            }
          }
        }

        return true;
      }


      bool _computeForceRaw(Real3D& force,
                            const Real3D& dist,
                            real distSqr) const {
        throw std::runtime_error("_computeForceRaw Not implemented!");
      }

      static LOG4ESPP_DECL_LOGGER(theLogger);
    };

    // provide pickle support
    struct LennardJonesLambda_pickle : boost::python::pickle_suite {
      static boost::python::tuple getinitargs(LennardJonesLambda const& pot) {
    	  real eps = pot.getEpsilon();
          real sig = pot.getSigma();
          real rc = pot.getCutoff();
          real sh = pot.getShift();
          real max_force = pot.getMaxForce();
          return boost::python::make_tuple(eps, sig, rc, sh, max_force);
      }
    };
  }
}

#endif
