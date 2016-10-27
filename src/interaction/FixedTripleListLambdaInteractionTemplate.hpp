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
#ifndef _INTERACTION_FIXEDTRIPLELISTLAMBDAINTERACTIONTEMPLATE_HPP
#define _INTERACTION_FIXEDTRIPLELISTLAMBDAINTERACTIONTEMPLATE_HPP

#include "mpi.hpp"
#include "Interaction.hpp"
#include "Real3D.hpp"
#include "Tensor.hpp"
#include "Particle.hpp"
#include "FixedTripleListLambda.hpp"
#include "esutil/Array2D.hpp"
#include "bc/BC.hpp"
#include "SystemAccess.hpp"
#include "types.hpp"

namespace espressopp {
namespace interaction {
template<typename _AngularPotential>
class FixedTripleListLambdaInteractionTemplate: public Interaction, SystemAccess {

 protected:
  typedef _AngularPotential Potential;

 public:
  FixedTripleListLambdaInteractionTemplate
      (shared_ptr<System> _system,
       shared_ptr<FixedTripleListLambda> _fixedtripleList,
       shared_ptr<Potential> _potential)
      : SystemAccess(_system), fixedtripleList(_fixedtripleList),
        potential(_potential) {
    if (!potential) {
      LOG4ESPP_ERROR(theLogger, "NULL potential");
    }
  }

  virtual ~FixedTripleListLambdaInteractionTemplate() {};

  void
  setFixedTripleList(shared_ptr<FixedTripleListLambda> _fixedtripleList) {
    fixedtripleList = _fixedtripleList;
  }

  shared_ptr<FixedTripleListLambda> getFixedTripleList() {
    return fixedtripleList;
  }

  void setPotential(shared_ptr<Potential> _potential) {
    if (_potential) {
      potential = _potential;
    } else {
      LOG4ESPP_ERROR(theLogger, "NULL potential");
    }
  }

  shared_ptr<Potential> getPotential() {
    return potential;
  }

  virtual void addForces();
  virtual real computeEnergy();
  virtual real computeEnergyDeriv();
  virtual real computeEnergyAA();
  virtual real computeEnergyCG();
  virtual void computeVirialX(std::vector<real> &p_xx_total, int bins);
  virtual real computeVirial();
  virtual void computeVirialTensor(Tensor &w);
  virtual void computeVirialTensor(Tensor &w, real z);
  virtual void computeVirialTensor(Tensor *w, int n);
  virtual real getMaxCutoff();
  virtual int bondType() { return Angular; }

 protected:
  int ntypes;
  shared_ptr<FixedTripleListLambda> fixedtripleList;
  //esutil::Array2D<Potential, esutil::enlarge> potentialArray;
  shared_ptr<Potential> potential;
};

//////////////////////////////////////////////////
// INLINE IMPLEMENTATION
//////////////////////////////////////////////////
template<typename _AngularPotential>
inline void
FixedTripleListLambdaInteractionTemplate<_AngularPotential>::
addForces() {
  LOG4ESPP_INFO(theLogger, "add forces computed by FixedTripleListLambda");
  const bc::BC &bc = *getSystemRef().bc;  // boundary conditions
  for (FixedTripleListLambda::Iterator it(*fixedtripleList); it.isValid(); ++it) {
    Particle &p1 = *it->first;
    Particle &p2 = *it->second;
    Particle &p3 = *it->third;
    real lambda = it->fourth;
    Real3D dist12, dist32;
    bc.getMinimumImageVectorBox(dist12, p1.position(), p2.position());
    bc.getMinimumImageVectorBox(dist32, p3.position(), p2.position());
    Real3D force12, force32;
    potential->_computeForce(force12, force32, dist12, dist32);
    p1.force() += lambda*force12;
    p2.force() -= lambda*force12 + lambda*force32;
    p3.force() += lambda*force32;
  }
}

template<typename _AngularPotential>
inline real
FixedTripleListLambdaInteractionTemplate<_AngularPotential>::
computeEnergy() {
  LOG4ESPP_INFO(theLogger, "compute energy of the triples");

  const bc::BC &bc = *getSystemRef().bc;
  real e = 0.0;
  for (FixedTripleListLambda::Iterator it(*fixedtripleList); it.isValid(); ++it) {
    const Particle &p1 = *it->first;
    const Particle &p2 = *it->second;
    const Particle &p3 = *it->third;
    const real lambda = it->fourth;
    Real3D dist12 = bc.getMinimumImageVector(p1.position(), p2.position());
    Real3D dist32 = bc.getMinimumImageVector(p3.position(), p2.position());
    e += lambda*potential->_computeEnergy(dist12, dist32);
  }
  real esum;
  boost::mpi::all_reduce(*mpiWorld, e, esum, std::plus<real>());
  return esum;
}

template<typename _AngularPotential>
inline real
FixedTripleListLambdaInteractionTemplate<_AngularPotential>::
computeEnergyDeriv() {
  std::cout << "Warning! At the moment computeEnergyDeriv() in FixedTripleListLambdaInteractionTemplate does not work."
            << std::endl;
  return 0.0;
}

template<typename _AngularPotential>
inline real
FixedTripleListLambdaInteractionTemplate<_AngularPotential>::
computeEnergyAA() {
  std::cout << "Warning! At the moment computeEnergyAA() in FixedTripleListLambdaInteractionTemplate does not work."
            << std::endl;
  return 0.0;
}

template<typename _AngularPotential>
inline real
FixedTripleListLambdaInteractionTemplate<_AngularPotential>::
computeEnergyCG() {
  std::cout << "Warning! At the moment computeEnergyCG() in FixedTripleListLambdaInteractionTemplate does not work."
            << std::endl;
  return 0.0;
}

template<typename _AngularPotential>
inline void
FixedTripleListLambdaInteractionTemplate<_AngularPotential>::
computeVirialX(std::vector<real> &p_xx_total, int bins) {
  std::cout << "Warning! At the moment computeVirialX in FixedTripleListLambdaInteractionTemplate does not work."
            << std::endl << "Therefore, the corresponding interactions won't be included in calculation." << std::endl;
}

template<typename _AngularPotential>
inline real
FixedTripleListLambdaInteractionTemplate<_AngularPotential>::
computeVirial() {
  LOG4ESPP_INFO(theLogger, "compute scalar virial of the triples");

  const bc::BC &bc = *getSystemRef().bc;
  real w = 0.0;
  for (FixedTripleListLambda::Iterator it(*fixedtripleList); it.isValid(); ++it) {
    const Particle &p1 = *it->first;
    const Particle &p2 = *it->second;
    const Particle &p3 = *it->third;
    const real lambda = it->fourth;
    //const Potential &potential = getPotential(p1.type(), p2.type());
    Real3D dist12, dist32;
    bc.getMinimumImageVectorBox(dist12, p1.position(), p2.position());
    bc.getMinimumImageVectorBox(dist32, p3.position(), p2.position());
    Real3D force12, force32;
    potential->_computeForce(force12, force32, dist12, dist32);
    w += dist12 * lambda*force12 + dist32 * lambda*force32;
  }
  real wsum;
  boost::mpi::all_reduce(*mpiWorld, w, wsum, std::plus<real>());
  return wsum;
}

template<typename _AngularPotential>
inline void
FixedTripleListLambdaInteractionTemplate<_AngularPotential>::
computeVirialTensor(Tensor &w) {
  LOG4ESPP_INFO(theLogger, "compute the virial tensor of the triples");

  Tensor wlocal(0.0);
  const bc::BC &bc = *getSystemRef().bc;
  for (FixedTripleListLambda::Iterator it(*fixedtripleList); it.isValid(); ++it) {
    const Particle &p1 = *it->first;
    const Particle &p2 = *it->second;
    const Particle &p3 = *it->third;
    const real lambda = it->fourth;
    //const Potential &potential = getPotential(0, 0);
    Real3D r12, r32;
    bc.getMinimumImageVectorBox(r12, p1.position(), p2.position());
    bc.getMinimumImageVectorBox(r32, p3.position(), p2.position());
    Real3D force12, force32;
    potential->_computeForce(force12, force32, r12, r32);
    wlocal += Tensor(r12, lambda*force12) + Tensor(r32, lambda*force32);
  }

  // reduce over all CPUs
  Tensor wsum(0.0);
  boost::mpi::all_reduce(*mpiWorld, (double *) &wlocal, 6, (double *) &wsum, std::plus<double>());
  w += wsum;
}

template<typename _AngularPotential>
inline void
FixedTripleListLambdaInteractionTemplate<_AngularPotential>::
computeVirialTensor(Tensor &w, real z) {
  LOG4ESPP_INFO(theLogger, "compute the virial tensor of the triples");

  std::cout << "Warning! At the moment IK computeVirialTensor for fixed triples does'n work" << std::endl;
}

template<typename _AngularPotential>
inline void
FixedTripleListLambdaInteractionTemplate<_AngularPotential>::
computeVirialTensor(Tensor *w, int n) {
  LOG4ESPP_INFO(theLogger, "compute the virial tensor of the triples");

  std::cout << "Warning! At the moment IK computeVirialTensor for fixed triples does'n work" << std::endl;
}

template<typename _AngularPotential>
inline real
FixedTripleListLambdaInteractionTemplate<_AngularPotential>::
getMaxCutoff() {
  return potential->getCutoff();
}
}
}
#endif
