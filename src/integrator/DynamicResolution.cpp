/*
  Copyright (C) 2015,2016
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

#include "DynamicResolution.hpp"
#include <algorithm>
#include <vector>
#include <iomanip>
#include <sstream>
#include "python.hpp"
#include "Real3D.hpp"
#include "Particle.hpp"
#include "Cell.hpp"
#include "System.hpp"
#include "storage/Storage.hpp"
#include "bc/BC.hpp"
#include "FixedTupleList.hpp"
#include "iterator/CellListAllPairsIterator.hpp"
#include "iterator/CellListIterator.hpp"

namespace espressopp {
namespace integrator {
using namespace espressopp::iterator;  //NOLINT

DynamicResolution::DynamicResolution(
    shared_ptr<System> _system,
    shared_ptr<TupleList> _fixedtupleList,
    real _rate) : Extension(_system),
                  fixedtupleList(_fixedtupleList) {
  LOG4ESPP_INFO(theLogger, "construct DynamicResolution");

  type = Extension::Adress;

  resolution_ = 0.0;
}

DynamicResolution::~DynamicResolution() {
  LOG4ESPP_INFO(theLogger, "~DynamicResolution");
  disconnect();
}

void DynamicResolution::connect() {
  _aftIntV = integrator->aftIntV.connect(
      boost::bind(&DynamicResolution::ChangeResolution, this),
      boost::signals2::at_back);
  // connection to after runInit()
  _SetPosVel = integrator->runInit.connect(
      boost::bind(&DynamicResolution::SetPosVel, this), boost::signals2::at_front);
  _befIntP = integrator->befIntP.connect(
      boost::bind(&DynamicResolution::SetPosVel, this), boost::signals2::at_front);
  // connection to inside of integrate1()
  _integrate1 = integrator->inIntP.connect(
      boost::bind(&DynamicResolution::integrate1, this, _1), boost::signals2::at_front);
  // connection to after _recalc2()
  _recalc2 = integrator->recalc2.connect(
      boost::bind(&DynamicResolution::aftCalcF, this), boost::signals2::at_front);
  // connection to after _befIntV()
  _befIntV = integrator->befIntV.connect(
      boost::bind(&DynamicResolution::aftCalcF, this), boost::signals2::at_front);
  // connection to after aftIntV()
  _aftIntV2 = integrator->aftIntV.connect(
      boost::bind(&DynamicResolution::SetVel, this), boost::signals2::at_front);
}

void DynamicResolution::disconnect() {
  _SetPosVel.disconnect();
  _recalc2.disconnect();
  _befIntV.disconnect();
  _integrate1.disconnect();
  _befIntP.disconnect();
  _aftIntV.disconnect();
  _aftIntV2.disconnect();
}

void DynamicResolution::set_active(bool active) {
  active_ = active;
}

void DynamicResolution::ChangeResolution() {
  if (!active_)
    return;

  // Increase the resolution lineary with the time.
  resolution_ += rate_;
  if (resolution_ > 1.0) {
    resolution_ = 1.0;
    active_ = false;
  } else if (resolution_ < 0.0) {
    resolution_ = 0.0;
    active_ = false;
  }
  updateWeights();
}


void DynamicResolution::integrate1(real &maxSqDist) {
  maxSqDist = SetPosVel();
}


void DynamicResolution::SetVel() {
  System &system = getSystemRef();
  // Set the positions and velocity of CG particles & update weights.
  real sqDist = 0.0;
  for (CellListIterator cit(system.storage->getLocalCells()); !cit.isDone(); ++cit) {
    Particle &vp = *cit;
    FixedTupleListAdress::iterator it3;
    it3 = fixedtupleList->find(&vp);

    if (it3 != fixedtupleList->end()) {
      std::vector<Particle *> atList;
      atList = it3->second;

      // Compute center of mass
      Real3D cmv(0.0, 0.0, 0.0);  // center of mass velocity

      for (std::vector<Particle *>::iterator it2 = atList.begin(); it2 != atList.end(); ++it2) {
        Particle &at = **it2;
        cmv += at.mass() * at.velocity();
        at.lambdaDeriv() = 0.0;
      }

      cmv /= vp.mass();

      // update (overwrite) the position and velocity of the VP
      vp.velocity() = cmv;

      vp.lambdaDeriv() = 0.0;
    }
  }
}

real DynamicResolution::SetPosVel() {
  System &system = getSystemRef();
  // Set the positions and velocity of CG particles & update weights.
  real sqDist = 0.0;
  for (CellListIterator cit(system.storage->getLocalCells()); !cit.isDone(); ++cit) {
    Particle &vp = *cit;
    FixedTupleListAdress::iterator it3;
    it3 = fixedtupleList->find(&vp);

    if (it3 != fixedtupleList->end()) {
      std::vector<Particle *> atList;
      atList = it3->second;

      // Compute center of mass
      Real3D cmp(0.0, 0.0, 0.0);  // center of mass position
      Real3D cmv(0.0, 0.0, 0.0);  // center of mass velocity

      for (std::vector<Particle *>::iterator it2 = atList.begin(); it2 != atList.end(); ++it2) {
        Particle &at = **it2;
        cmp += at.mass() * at.position();
        cmv += at.mass() * at.velocity();
        at.lambdaDeriv() = 0.0;
        at.lambda() = resolution_;
      }

      cmp /= vp.mass();
      cmv /= vp.mass();

      // Fold position.
      system.bc->foldPosition(cmp);

      // update (overwrite) the position and velocity of the VP
      Real3D deltaP = vp.position() - cmp;
      sqDist = std::max(sqDist, deltaP.sqr());
      vp.position() = cmp;
      vp.velocity() = cmv;

      vp.lambda() = resolution_;
      vp.lambdaDeriv() = 0.0;
    }
  }
  return sqDist;
}

void DynamicResolution::updateWeights() {
  System &system = getSystemRef();
  for (CellListIterator cit(system.storage->getLocalCells()); !cit.isDone(); ++cit) {
    Particle &vp = *cit;
    FixedTupleListAdress::iterator it3;
    it3 = fixedtupleList->find(&vp);

    vp.lambda() = resolution_;
    vp.lambdaDeriv() = 0.0;

    if (it3 != fixedtupleList->end()) {
      std::vector<Particle *> atList;
      atList = it3->second;

      // Propagate lambda/lambdaDeriv downstream to underlying atoms
      for (std::vector<Particle *>::iterator it2 = atList.begin(); it2 != atList.end(); ++it2) {
        Particle &at = **it2;
        at.lambda() = resolution_;
        at.lambdaDeriv() = 0.0;
      }
    }
  }
}

void DynamicResolution::aftCalcF() {
  System &system = getSystemRef();
  for (CellListIterator cit(system.storage->getLocalCells()); !cit.isDone(); ++cit) {
    Particle &vp = *cit;

    FixedTupleListAdress::iterator it3;
    it3 = fixedtupleList->find(&vp);

    if (it3 != fixedtupleList->end()) {
      std::vector<Particle *> atList;
      atList = it3->second;

      // update force of AT particles belonging to a VP
      Real3D vpfm = vp.force() / vp.getMass();
      for (std::vector<Particle *>::iterator it2 = atList.begin(); it2 != atList.end(); ++it2) {
        Particle &at = **it2;

        at.force() += at.mass() * vpfm;
      }
    }
  }
}

void DynamicResolution::registerPython() {
  using namespace espressopp::python;  // NOLINT
  class_<DynamicResolution, shared_ptr<DynamicResolution>, bases<Extension> >
      ("integrator_DynamicResolution", init<shared_ptr<System>,
                                            shared_ptr<FixedTupleList>, real>())
      .add_property("active", &DynamicResolution::active, &DynamicResolution::set_active)
      .add_property("resolution",
                    &DynamicResolution::resolution,
                    &DynamicResolution::set_resolution)
      .add_property("rate", &DynamicResolution::rate, &DynamicResolution::set_rate)
      .def("connect", &DynamicResolution::connect)
      .def("disconnect", &DynamicResolution::disconnect)
      .def("SetPosVel", &DynamicResolution::SetPosVel);
}


BasicDynamicResolutionType::BasicDynamicResolutionType(shared_ptr<System> system):
    Extension(system) {
  type = Extension::Adress;
}

BasicDynamicResolutionType::~BasicDynamicResolutionType() {
  LOG4ESPP_INFO(theLogger, "~DynamicResolution");
  disconnect();
}

void BasicDynamicResolutionType::connect() {
  _aftIntV = integrator->aftIntV.connect(
      boost::bind(&BasicDynamicResolutionType::UpdateWeights, this),
      boost::signals2::at_back);
}

void BasicDynamicResolutionType::disconnect() {
  _aftIntV.disconnect();
}

void BasicDynamicResolutionType::registerPython() {
  using namespace espressopp::python;  // NOLINT
  class_<BasicDynamicResolutionType, shared_ptr<BasicDynamicResolutionType>, bases<Extension> >
      ("integrator_BasicDynamicResolutionType", init<shared_ptr<System> >())
      .def("set_type_rate", &BasicDynamicResolutionType::setTypeRate)
       .def("connect", &BasicDynamicResolutionType::connect)
       .def("disconnect", &BasicDynamicResolutionType::disconnect)
       .def("add_postprocess", &BasicDynamicResolutionType::addPostProcess);
}

void BasicDynamicResolutionType::UpdateWeights() {
  real time0 = wallTimer.getElapsedTime();
  System &system = getSystemRef();
  for (CellListIterator cit(system.storage->getLocalCells()); !cit.isDone(); ++cit) {
    Particle &vp = *cit;

    real new_lambda = vp.lambda() + rate_type_[vp.type()];
    bool lambda_0(false);
    bool lambda_1(false);
    if (new_lambda < 0.0) {
      new_lambda = 0.0;
      lambda_0 = true;
    } else if (new_lambda > 1.0) {
      new_lambda = 1.0;
      lambda_1 = true;
    }

    vp.lambda() = new_lambda;
    vp.lambdaDeriv() = 0.0;

    /** Run PostProcess methods whenever resolution reaches 1 or 0. */
    if (lambda_0) {
      for (PostProcessMap::iterator it = post_process_0.begin(); it != post_process_0.end(); ++it) {
        (*it)->process(vp, vp);
      }
    } else if (lambda_1) {
      for (PostProcessMap::iterator it = post_process_1.begin(); it != post_process_1.end(); ++it) {
        (*it)->process(vp, vp);
      }
    }
  }
  timeUpdateWeights += wallTimer.getElapsedTime() - time0;
}

void FixedListDynamicResolution::connect() {
  _aftIntV = integrator->aftIntV.connect(boost::bind(&FixedListDynamicResolution::updateLists, this));
}

void FixedListDynamicResolution::disconnect() {
  _aftIntV.disconnect();
}

/**
 * Update lambda parameters on fixed lists.
 */
void FixedListDynamicResolution::updateLists() {
  real time0 = wallTimer.getElapsedTime();
  // Update FixedPairLists
  for (FixedPairListRate::iterator it = fixed_pair_list_rate_.begin(); it != fixed_pair_list_rate_.end(); ++it) {
    it->first->incrementAllLambda(it->second);
  }

  // Update TripleLists
  for (FixedTripleListRate::iterator it = fixed_triple_list_rate_.begin(); it != fixed_triple_list_rate_.end(); ++it) {
    it->first->incrementAllLambda(it->second);
  }

  // Update QuadrupleLists
  for (FixedQuadrupleListRate::iterator it = fixed_quadruple_list_rate_.begin();
       it != fixed_quadruple_list_rate_.end(); ++it) {
    it->first->incrementAllLambda(it->second);
  }
  timeUpdateLists += wallTimer.getElapsedTime() - time0;
}

void FixedListDynamicResolution::registerPython() {
  using namespace espressopp::python;  // NOLINT
  class_<FixedListDynamicResolution, shared_ptr<FixedListDynamicResolution>, bases<Extension> >
      ("integrator_FixedListDynamicResolution", init<shared_ptr<System> >())
      .def("connect", &FixedListDynamicResolution::connect)
      .def("disconnect", &FixedListDynamicResolution::disconnect)
      .def("update_lists", &FixedListDynamicResolution::updateLists)
      .def("register_pair_list", &FixedListDynamicResolution::registerPairList)
      .def("register_triple_list", &FixedListDynamicResolution::registerTripleList)
      .def("register_quadruple_list", &FixedListDynamicResolution::registerQuadrupleList)
      .def("update_pair_list", &FixedListDynamicResolution::updatePairList)
      .def("update_triple_list", &FixedListDynamicResolution::updateTripleList)
      .def("update_quadruple_list", &FixedListDynamicResolution::updateQuadrupleList);
}

}  // end namespace integrator
}  // end namespace espressopp
