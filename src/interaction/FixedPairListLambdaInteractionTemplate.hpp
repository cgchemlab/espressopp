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
#ifndef _INTERACTION_FIXEDPAIRLISTLAMBDAINTERACTIONTEMPLATE_HPP
#define _INTERACTION_FIXEDPAIRLISTLAMBDAINTERACTIONTEMPLATE_HPP

#include "mpi.hpp"
#include "Interaction.hpp"
#include "Real3D.hpp"
#include "Tensor.hpp"
#include "Particle.hpp"
#include "FixedPairListLambda.hpp"
#include "esutil/Array2D.hpp"
#include "bc/BC.hpp"
#include "SystemAccess.hpp"
#include "Interaction.hpp"
#include "types.hpp"

namespace espressopp {
namespace interaction {
template<typename _Potential>
class FixedPairListLambdaInteractionTemplate: public Interaction, SystemAccess {

 protected:
  typedef _Potential Potential;

 public:
  FixedPairListLambdaInteractionTemplate
      (shared_ptr <System> system,
       shared_ptr <FixedPairListLambda> _fixedpairList,
       shared_ptr <Potential> _potential)
      : SystemAccess(system), fixedpairList(_fixedpairList),
        potential(_potential) {
    if (!potential) {
      LOG4ESPP_ERROR(theLogger, "NULL potential");
    }
  }

  virtual ~FixedPairListLambdaInteractionTemplate() {};

  void
  setFixedPairList(shared_ptr <FixedPairListLambda> _fixedpairList) {
    fixedpairList = _fixedpairList;
  }

  shared_ptr <FixedPairListLambda> getFixedPairList() {
    return fixedpairList;
  }

  void
  setPotential(shared_ptr <Potential> _potential) {
    if (_potential) {
      potential = _potential;
    } else {
      LOG4ESPP_ERROR(theLogger, "NULL potential");
    }
  }

  shared_ptr <Potential> getPotential() {
    return potential;
  }

  virtual void addForces();
  virtual real computeEnergy();
  virtual real computeEnergyDeriv();
  virtual real computeEnergyAA();
  virtual real computeEnergyCG();
  virtual void computeVirialX(std::vector <real> &p_xx_total, int bins);
  virtual real computeVirial();
  virtual void computeVirialTensor(Tensor &w);
  virtual void computeVirialTensor(Tensor &w, real z);
  virtual void computeVirialTensor(Tensor *w, int n);
  virtual real getMaxCutoff();
  virtual int bondType() { return Pair; }

 protected:
  int ntypes;
  shared_ptr <FixedPairListLambda> fixedpairList;
  shared_ptr <Potential> potential;
};

//////////////////////////////////////////////////
// INLINE IMPLEMENTATION
//////////////////////////////////////////////////
template<typename _Potential>
inline void
FixedPairListLambdaInteractionTemplate<_Potential>::addForces() {
  LOG4ESPP_INFO(_Potential::theLogger, "adding forces of FixedPairListLambda");
  const bc::BC &bc = *getSystemRef().bc;  // boundary conditions
  real ltMaxBondSqr = fixedpairList->getLongtimeMaxBondSqr();
  for (FixedPairListLambda::Iterator it(*fixedpairList); it.isValid(); ++it) {
    Particle &p1 = *it->first;
    Particle &p2 = *it->second;
    real lambda = it->third;
    Real3D dist;
    bc.getMinimumImageVectorBox(dist, p1.position(), p2.position());
    Real3D force;
    real d = dist.sqr();
    if (d > ltMaxBondSqr) {
      fixedpairList->setLongtimeMaxBondSqr(d);
      ltMaxBondSqr = d;
    }
    if (potential->_computeForce(force, dist)) {
      p1.force() += lambda*force;
      p2.force() -= lambda*force;
      LOG4ESPP_DEBUG(_Potential::theLogger,
                     "p" << p1.id() << "(" << p1.position()[0] << "," << p1.position()[1] << "," << p1.position()[2]
                         << ") "
                         << "p" << p2.id() << "(" << p2.position()[0] << "," << p2.position()[1] << ","
                         << p2.position()[2] << ") "
                         << "dist=" << sqrt(dist * dist) << " "
                         << "force=(" << force[0] << "," << force[1] << "," << force[2] << ")"
                         << " lambda=" << lambda
                    );
    }
  }
}

template<typename _Potential>
inline real FixedPairListLambdaInteractionTemplate<_Potential>::computeEnergy() {

  LOG4ESPP_INFO(theLogger, "compute energy of the FixedPairListLambda pairs");

  real e = 0.0;
  const bc::BC &bc = *getSystemRef().bc;  // boundary conditions
  for (FixedPairListLambda::Iterator it(*fixedpairList);
       it.isValid(); ++it) {
    const Particle &p1 = *it->first;
    const Particle &p2 = *it->second;
    real lambda = it->third;
    Real3D r21;
    bc.getMinimumImageVectorBox(r21, p1.position(), p2.position());
    e += lambda*potential->_computeEnergy(r21);
  }
  real esum;
  boost::mpi::all_reduce(*mpiWorld, e, esum, std::plus<real>());
  return esum;
}

template<typename _Potential>
inline real FixedPairListLambdaInteractionTemplate<_Potential>::computeEnergyDeriv() {
  std::cout << "Warning! At the moment computeEnergyDeriv() in FixedPairListLambdaInteractionTemplate does not work."
            << std::endl;
  return 0.0;
}

template<typename _Potential>
inline real FixedPairListLambdaInteractionTemplate<_Potential>::computeEnergyAA() {
  std::cout << "Warning! At the moment computeEnergyAA() in FixedPairListLambdaInteractionTemplate does not work."
            << std::endl;
  return 0.0;
}

template<typename _Potential>
inline real FixedPairListLambdaInteractionTemplate<_Potential>::computeEnergyCG() {
  std::cout << "Warning! At the moment computeEnergyCG() in FixedPairListLambdaInteractionTemplate does not work."
            << std::endl;
  return 0.0;
}

template<typename _Potential>
inline void FixedPairListLambdaInteractionTemplate<_Potential>::computeVirialX(std::vector <real> &p_xx_total, int bins) {
  LOG4ESPP_INFO(theLogger, "compute virial p_xx of the pressure tensor slabwise");
}


template<typename _Potential>
inline real FixedPairListLambdaInteractionTemplate<_Potential>::computeVirial() {
  LOG4ESPP_INFO(theLogger, "compute the virial for the FixedPair List");

  real w = 0.0;
  const bc::BC &bc = *getSystemRef().bc;  // boundary conditions
  for (FixedPairListLambda::Iterator it(*fixedpairList);
       it.isValid(); ++it) {
    const Particle &p1 = *it->first;
    const Particle &p2 = *it->second;
    real lambda = it->third;

    Real3D r21;
    bc.getMinimumImageVectorBox(r21, p1.position(), p2.position());
    Real3D force;
    if (potential->_computeForce(force, r21)) {
      w += r21 * (lambda*force);
    }
  }

  real wsum;
  boost::mpi::all_reduce(*mpiWorld, w, wsum, std::plus<real>());
  return wsum;
}

template<typename _Potential>
inline void FixedPairListLambdaInteractionTemplate<_Potential>::computeVirialTensor(Tensor &w) {
  LOG4ESPP_INFO(theLogger, "compute the virial tensor for the FixedPair List");

  Tensor wlocal(0.0);
  const bc::BC &bc = *getSystemRef().bc;  // boundary conditions
  for (FixedPairListLambda::Iterator it(*fixedpairList);
       it.isValid(); ++it) {
    const Particle &p1 = *it->first;
    const Particle &p2 = *it->second;
    real lambda = it->third;

    Real3D r21;
    bc.getMinimumImageVectorBox(r21, p1.position(), p2.position());
    Real3D force;
    if (potential->_computeForce(force, r21)) {
      wlocal += Tensor(r21, lambda*force);
    }
  }

  // reduce over all CPUs
  Tensor wsum(0.0);
  boost::mpi::all_reduce(*mpiWorld, (double *) &wlocal, 6, (double *) &wsum, std::plus<double>());
  w += wsum;
}

template<typename _Potential>
inline void FixedPairListLambdaInteractionTemplate<_Potential>::computeVirialTensor(Tensor &w, real z) {

}

template<typename _Potential>
inline void FixedPairListLambdaInteractionTemplate<_Potential>::computeVirialTensor(Tensor *w, int n) {

}

template<typename _Potential>
inline real FixedPairListLambdaInteractionTemplate<_Potential>::
getMaxCutoff() {
  return potential->getCutoff();
}
}
}
#endif
