/*
  Copyright (c) 2015
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
#ifndef _DYNAMICRESOLUTION_HPP
#define _DYNAMICRESOLUTION_HPP

#include "log4espp.hpp"
#include "types.hpp"
#include "Particle.hpp"
#include "SystemAccess.hpp"
#include "VerletListAdress.hpp"
#include "FixedTupleListAdress.hpp"
#include "Extension.hpp"
#include "VelocityVerlet.hpp"
#include "Adress.hpp"

#include "boost/signals2.hpp"


namespace espressopp {
namespace integrator {

/**
 * Simple extension that change only lambda parameter (0.0 - 1.0)
 * depends on the particle type.
 */
class BasicDynamicResolutionType : public Extension {
 public:
  BasicDynamicResolutionType(shared_ptr<System> _system);
  ~BasicDynamicResolutionType();

  void SetTypeRate(longint type, real rate) { rate_type_[type] = rate; }

  static void registerPython();
 private:
  void connect();
  void disconnect();
  boost::unordered_map<longint, real> rate_type_;
  void UpdateWeights();
  boost::signals2::connection _aftIntV;
};

/**
 * This module implement dynamic resolution extension.
 */
class DynamicResolution : public Extension {
 public:
  DynamicResolution(shared_ptr<System> _system,
    shared_ptr<TupleList> _fixedtupleList,
    real _rate);

  ~DynamicResolution();

  real resolution() { return resolution_; }
  void set_resolution(real resolution) { resolution_ = resolution;}

  real rate() { return rate_; }
  void set_rate(real val) { rate_ = val; }

  bool active() { return active_; }
  void set_active(bool active);

  /** Register this class so it can be used from Python. */
  static void registerPython();

  /** Static method for compute weights. */
  static real ComputeWeight(real w1, real w2) {
    return w1*w1;
  }

  static real ComputeWeight(real w1, real w2, real w3) {
    return pow(w1*w2*w3, 2.0/3.0);
  }

  static real ComputeWeight(real w1, real w2, real w3, real w4) {
    return pow(w1*w2*w3*w4, 1.0/2.0);
  }

 private:
  real weight(real) { return resolution_; }
  real weightderivative(real) { return 0.0; }
  void updateWeights();
  real SetPosVel();
  void SetVel();
  void integrate1(real &maxSqDist);
  void aftCalcF();

  void ChangeResolution();
  real rate_;
  bool active_;

  shared_ptr<TupleList> fixedtupleList;

  real resolution_;  /// Current value of resolution, between 0.0 - 1.0
  void connect();
  void disconnect();

  // Signals
  boost::signals2::connection _aftIntV, _SetPosVel, _befIntP, _integrate1;
  boost::signals2::connection _recalc2, _befIntV, _aftIntV2;
};
}  // namespace integrator
}  // namespace espressopp
#endif
