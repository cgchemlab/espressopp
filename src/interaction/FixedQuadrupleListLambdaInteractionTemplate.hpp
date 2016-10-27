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
#ifndef _INTERACTION_FIXEDQUADRUPLELISTLAMBDAINTERACTIONTEMPLATE_HPP
#define _INTERACTION_FIXEDQUADRUPLELISTLAMBDAINTERACTIONTEMPLATE_HPP

#include "mpi.hpp"
#include "Interaction.hpp"
#include "Real3D.hpp"
#include "Tensor.hpp"
#include "Particle.hpp"
#include "FixedQuadrupleListLambda.hpp"
#include "esutil/Array2D.hpp"
#include "bc/BC.hpp"
#include "SystemAccess.hpp"
#include "types.hpp"

namespace espressopp {
namespace interaction {
template<typename _DihedralPotential>
class FixedQuadrupleListLambdaInteractionTemplate: public Interaction, SystemAccess {
 protected:
  typedef _DihedralPotential Potential;
 public:
  FixedQuadrupleListLambdaInteractionTemplate
      (shared_ptr<System> _system,
       shared_ptr<FixedQuadrupleListLambda> _fixedquadrupleList,
       shared_ptr<Potential> _potential)
      : SystemAccess(_system), fixedquadrupleList(_fixedquadrupleList),
        potential(_potential) {
    if (!potential) {
      LOG4ESPP_ERROR(theLogger, "NULL potential");
    }

  }

  void
  setFixedQuadrupleList(shared_ptr<FixedQuadrupleListLambda> _fixedquadrupleList) {
    fixedquadrupleList = _fixedquadrupleList;
  }

  virtual ~FixedQuadrupleListLambdaInteractionTemplate() {};

  shared_ptr<FixedQuadrupleListLambda> getFixedQuadrupleList() {
    return fixedquadrupleList;
  }

  void
  setPotential(shared_ptr<Potential> _potential) {
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
  virtual int bondType() { return Dihedral; }

 protected:
  int ntypes;
  shared_ptr<FixedQuadrupleListLambda> fixedquadrupleList;
  shared_ptr<Potential> potential;
};

//////////////////////////////////////////////////
// INLINE IMPLEMENTATION
//////////////////////////////////////////////////
template<typename _DihedralPotential>
inline void
FixedQuadrupleListLambdaInteractionTemplate<_DihedralPotential>::
addForces() {

  LOG4ESPP_INFO(theLogger, "add forces computed by FixedQuadrupleList");

  const bc::BC &bc = *getSystemRef().bc;  // boundary conditions

  for (FixedQuadrupleListLambda::Iterator it(*fixedquadrupleList); it.isValid(); ++it) {
    Particle &p1 = *it->first;
    Particle &p2 = *it->second;
    Particle &p3 = *it->third;
    Particle &p4 = *it->fourth;
    real lambda = it->fifth;

    Real3D dist21, dist32, dist43; //

    bc.getMinimumImageVectorBox(dist21, p2.position(), p1.position());
    bc.getMinimumImageVectorBox(dist32, p3.position(), p2.position());
    bc.getMinimumImageVectorBox(dist43, p4.position(), p3.position());

    Real3D force1, force2, force3, force4;  // result forces

    potential->_computeForce(force1, force2, force3, force4,
                             dist21, dist32, dist43);
    p1.force() += lambda*force1;
    p2.force() += lambda*force2;
    p3.force() += lambda*force3;
    p4.force() += lambda*force4;
  }
}

template<typename _DihedralPotential>
inline real
FixedQuadrupleListLambdaInteractionTemplate<_DihedralPotential>::
computeEnergy() {
  LOG4ESPP_INFO(theLogger, "compute energy of the quadruples");

  const bc::BC &bc = *getSystemRef().bc;  // boundary conditions
  real e = 0.0;
  for (FixedQuadrupleListLambda::Iterator it(*fixedquadrupleList); it.isValid(); ++it) {
    const Particle &p1 = *it->first;
    const Particle &p2 = *it->second;
    const Particle &p3 = *it->third;
    const Particle &p4 = *it->fourth;
    real lambda = it->fifth;

    Real3D dist21, dist32, dist43; //

    bc.getMinimumImageVectorBox(dist21, p2.position(), p1.position());
    bc.getMinimumImageVectorBox(dist32, p3.position(), p2.position());
    bc.getMinimumImageVectorBox(dist43, p4.position(), p3.position());

    e += lambda*potential->_computeEnergy(dist21, dist32, dist43);
  }
  real esum;
  boost::mpi::all_reduce(*mpiWorld, e, esum, std::plus<real>());
  return esum;
}

template<typename _DihedralPotential>
inline real
FixedQuadrupleListLambdaInteractionTemplate<_DihedralPotential>::
computeEnergyDeriv() {
  std::cout
      << "Warning! At the moment computeEnergyDeriv() in FixedQuadrupleListLambdaInteractionTemplate does not work."
      << std::endl;
  return 0.0;
}

template<typename _DihedralPotential>
inline real
FixedQuadrupleListLambdaInteractionTemplate<_DihedralPotential>::
computeEnergyAA() {
  std::cout << "Warning! At the moment computeEnergyAA() in FixedQuadrupleListLambdaInteractionTemplate does not work."
            << std::endl;
  return 0.0;
}

template<typename _DihedralPotential>
inline real
FixedQuadrupleListLambdaInteractionTemplate<_DihedralPotential>::
computeEnergyCG() {
  std::cout << "Warning! At the moment computeEnergyCG() in FixedQuadrupleListLambdaInteractionTemplate does not work."
            << std::endl;
  return 0.0;
}

template<typename _DihedralPotential>
inline void
FixedQuadrupleListLambdaInteractionTemplate<_DihedralPotential>::
computeVirialX(std::vector<real> &p_xx_total, int bins) {
  std::cout << "Warning! At the moment computeVirialX in FixedQuadrupleListLambdaInteractionTemplate does not work."
            << std::endl << "Therefore, the corresponding interactions won't be included in calculation." << std::endl;
}

template<typename _DihedralPotential>
inline real
FixedQuadrupleListLambdaInteractionTemplate<_DihedralPotential>::
computeVirial() {
  LOG4ESPP_INFO(theLogger, "compute scalar virial of the quadruples");

  real w = 0.0;
  const bc::BC &bc = *getSystemRef().bc;  // boundary conditions
  for (FixedQuadrupleListLambda::Iterator it(*fixedquadrupleList); it.isValid(); ++it) {
    const Particle &p1 = *it->first;
    const Particle &p2 = *it->second;
    const Particle &p3 = *it->third;
    const Particle &p4 = *it->fourth;
    real lambda = it->fifth;

    Real3D dist21, dist32, dist43;

    bc.getMinimumImageVectorBox(dist21, p2.position(), p1.position());
    bc.getMinimumImageVectorBox(dist32, p3.position(), p2.position());
    bc.getMinimumImageVectorBox(dist43, p4.position(), p3.position());

    Real3D force1, force2, force3, force4;

    potential->_computeForce(force1, force2, force3, force4,
                             dist21, dist32, dist43);

    w += dist21 * lambda*force1 + dist32 * lambda*force2;
  }

  real wsum;
  boost::mpi::all_reduce(*mpiWorld, w, wsum, std::plus<real>());
  return w;
}

template<typename _DihedralPotential>
inline void
FixedQuadrupleListLambdaInteractionTemplate<_DihedralPotential>::
computeVirialTensor(Tensor &w) {
  LOG4ESPP_INFO(theLogger, "compute the virial tensor of the quadruples");

  Tensor wlocal(0.0);
  const bc::BC &bc = *getSystemRef().bc;

  for (FixedQuadrupleListLambda::Iterator it(*fixedquadrupleList); it.isValid(); ++it) {
    const Particle &p1 = *it->first;
    const Particle &p2 = *it->second;
    const Particle &p3 = *it->third;
    const Particle &p4 = *it->fourth;
    real lambda = it->fifth;

    Real3D dist21, dist32, dist43;

    bc.getMinimumImageVectorBox(dist21, p2.position(), p1.position());
    bc.getMinimumImageVectorBox(dist32, p3.position(), p2.position());
    bc.getMinimumImageVectorBox(dist43, p4.position(), p3.position());

    Real3D force1, force2, force3, force4;

    potential->_computeForce(force1, force2, force3, force4,
                             dist21, dist32, dist43);

    wlocal += Tensor(dist21, lambda*force1) - Tensor(dist32, lambda*force2);
  }
  // reduce over all CPUs
  Tensor wsum(0.0);
  boost::mpi::all_reduce(*mpiWorld, (double *) &wlocal, 6, (double *) &wsum, std::plus<double>());
  w += wsum;
}


template<typename _DihedralPotential>
inline void
FixedQuadrupleListLambdaInteractionTemplate<_DihedralPotential>::
computeVirialTensor(Tensor &w, real z) {
  LOG4ESPP_INFO(theLogger, "compute the virial tensor of the quadruples");

  Tensor wlocal(0.0);
  const bc::BC &bc = *getSystemRef().bc;

  std::cout << "Warning!!! computeVirialTensor in specified volume doesn't work for "
      "FixedQuadrupleListLambdaLambdaInteractionTemplate at the moment" << std::endl;

  for (FixedQuadrupleListLambda::Iterator it(*fixedquadrupleList); it.isValid(); ++it) {
    const Particle &p1 = *it->first;
    const Particle &p2 = *it->second;
    const Particle &p3 = *it->third;
    const Particle &p4 = *it->fourth;
    real lambda = it->fifth;

    Real3D dist21, dist32, dist43;

    bc.getMinimumImageVectorBox(dist21, p2.position(), p1.position());
    bc.getMinimumImageVectorBox(dist32, p3.position(), p2.position());
    bc.getMinimumImageVectorBox(dist43, p4.position(), p3.position());

    Real3D force1, force2, force3, force4;

    potential->_computeForce(force1, force2, force3, force4,
                             dist21, dist32, dist43);

    wlocal += Tensor(dist21, lambda*force1) - Tensor(dist32, lambda*force2);
  }
  // reduce over all CPUs
  Tensor wsum(0.0);
  boost::mpi::all_reduce(*mpiWorld, (double *) &wlocal, 6, (double *) &wsum, std::plus<double>());
  w += wsum;
}

template<typename _DihedralPotential>
inline void
FixedQuadrupleListLambdaInteractionTemplate<_DihedralPotential>::
computeVirialTensor(Tensor *w, int n) {
  std::cout << "Warning!!! computeVirialTensor in specified volume doesn't work for "
      "FixedQuadrupleListLambdaInteractionTemplate at the moment" << std::endl;
}

template<typename _DihedralPotential>
inline real
FixedQuadrupleListLambdaInteractionTemplate<_DihedralPotential>::
getMaxCutoff() {
  return potential->getCutoff();
}
}
}
#endif
