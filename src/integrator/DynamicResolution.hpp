/*
  Copyright (c) 2015,2016
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
#include "ChemicalReactionPostProcess.hpp"
#include "FixedPairListLambda.hpp"
#include "FixedTripleListLambda.hpp"
#include "FixedQuadrupleListLambda.hpp"

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

  void setTypeRate(longint type, real rate) { rate_type_[type] = rate; }
  /**
   * Define action after lambda reaches 1.0 or 0.0.
   */
  void addPostProcess(const shared_ptr<integrator::ChemicalReactionPostProcess> pp, int when = 1) {
    if (when == 1) {
      post_process_1.push_back(pp);
    } else if (when == 0) {
      post_process_0.push_back(pp);
    } else {
      throw std::runtime_error("Wrong value of when.");
    }
  }

  static void registerPython();
 private:
  void connect();
  void disconnect();
  boost::unordered_map<longint, real> rate_type_;
  void UpdateWeights();
  boost::signals2::connection _aftIntV;
  std::vector<shared_ptr<integrator::ChemicalReactionPostProcess> > post_process_0;
  std::vector<shared_ptr<integrator::ChemicalReactionPostProcess> > post_process_1;
};


class FixedListDynamicResolution : public Extension {
 public:
  FixedListDynamicResolution(shared_ptr<System> system_)
      : Extension(system_) {}
  ~FixedListDynamicResolution() {
    disconnect();
  }

  void registerPairList(shared_ptr<FixedPairListLambda> fpl, real rate) {
    fixed_pair_list_rate_.push_back(std::make_pair(fpl, rate));
  }
  void registerTripleList(shared_ptr<FixedTripleListLambda> ftl, real rate) {
    fixed_triple_list_rate_.push_back(std::make_pair(ftl, rate));
  }
  void registerQuadrupleList(shared_ptr<FixedQuadrupleListLambda> fql, real rate) {
    fixed_quadruple_list_rate_.push_back(std::make_pair(fql, rate));
  }

  void updateLists();

  static void registerPython();
 private:
  void connect();
  void disconnect();

  typedef std::vector<std::pair<shared_ptr<FixedPairListLambda>, real> > FixedPairListRate;
  typedef std::vector<std::pair<shared_ptr<FixedTripleListLambda>, real> > FixedTripleListRate;
  typedef std::vector<std::pair<shared_ptr<FixedQuadrupleListLambda>, real> > FixedQuadrupleListRate;

  FixedPairListRate fixed_pair_list_rate_;
  FixedTripleListRate fixed_triple_list_rate_;
  FixedQuadrupleListRate fixed_quadruple_list_rate_;

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
