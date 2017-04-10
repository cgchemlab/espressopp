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
#ifndef _INTERACTION_FIXEDQUADRUPLELISTTYPESLAMBDAINTERACTIONTEMPLATE_HPP
#define _INTERACTION_FIXEDQUADRUPLELISTTYPESLAMBDAINTERACTIONTEMPLATE_HPP

#include <algorithm>
#include <functional>
#include <vector>

#include "mpi.hpp"
#include "Interaction.hpp"
#include "Real3D.hpp"
#include "Tensor.hpp"
#include "Particle.hpp"
#include "FixedQuadrupleListLambda.hpp"
#include "esutil/Array4D.hpp"
#include "bc/BC.hpp"
#include "SystemAccess.hpp"
#include "types.hpp"

namespace espressopp {
namespace interaction {
template<typename _DihedralPotential>
class FixedQuadrupleListTypesLambdaInteractionTemplate: public Interaction, SystemAccess {
 protected:
  typedef _DihedralPotential Potential;

 public:
  FixedQuadrupleListTypesLambdaInteractionTemplate
      (shared_ptr<System> _system,
       shared_ptr<FixedQuadrupleListLambda> _fixedquadrupleList)
      : SystemAccess(_system), fixedquadrupleList(_fixedquadrupleList) {
    potentialArray = esutil::Array4D<Potential, esutil::enlarge>(0, 0, 0, 0, Potential());
    ntypes = 0;
  }

  void
  setFixedQuadrupleList(shared_ptr<FixedQuadrupleListLambda> _fixedquadrupleList) {
    fixedquadrupleList = _fixedquadrupleList;
  }

  virtual ~FixedQuadrupleListTypesLambdaInteractionTemplate() { }

  shared_ptr<FixedQuadrupleListLambda> getFixedQuadrupleList() {
    return fixedquadrupleList;
  }

  void setPotential(int type1, int type2, int type3, int type4, const Potential &potential) {
    // typeX+1 because i<ntypes
    ntypes = std::max(ntypes, std::max(std::max(std::max(type1 + 1, type2 + 1), type3+1), type4+1));
    potentialArray.at(type1, type2, type3, type4) = potential;
    potentialArray.at(type4, type3, type2, type1) = potential;
  }

  // this is used in the innermost force-loop
  Potential &getPotential(int type1, int type2, int type3, int type4) {
    return potentialArray.at(type1, type2, type3, type4);
  }

  shared_ptr<Potential> getPotentialPtr(int type1, int type2, int type3, int type4) {
    return make_shared<Potential>(potentialArray.at(type1, type2, type3, type4));
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
  esutil::Array4D<Potential, esutil::enlarge> potentialArray;
};

//////////////////////////////////////////////////
// INLINE IMPLEMENTATION
//////////////////////////////////////////////////
template<typename _DihedralPotential>
inline void
FixedQuadrupleListTypesLambdaInteractionTemplate<_DihedralPotential>::
addForces() {
  LOG4ESPP_INFO(theLogger, "add forces computed by FixedQuadrupleList");

  const bc::BC &bc = *getSystemRef().bc;  // boundary conditions

  for (FixedQuadrupleListLambda::IteratorParticleLambda it(fixedquadrupleList->getParticleQuadruples()); it.isValid(); ++it) {
    Particle &p1 = *it->first;
    Particle &p2 = *it->second;
    Particle &p3 = *it->third;
    Particle &p4 = *it->fourth;
    real lambda = it->lambda;

    longint type1 = p1.type();
    longint type2 = p2.type();
    longint type3 = p3.type();
    longint type4 = p4.type();

    const Potential &potential = getPotential(type1, type2, type3, type4);

    Real3D dist21, dist32, dist43;

    bc.getMinimumImageVectorBox(dist21, p2.position(), p1.position());
    bc.getMinimumImageVectorBox(dist32, p3.position(), p2.position());
    bc.getMinimumImageVectorBox(dist43, p4.position(), p3.position());

    Real3D force1, force2, force3, force4;  // result forces

    potential.computeForce(force1, force2, force3, force4,
                           dist21, dist32, dist43);
    p1.force() += lambda*force1;
    p2.force() += lambda*force2;
    p3.force() += lambda*force3;
    p4.force() += lambda*force4;
  }
}

template<typename _DihedralPotential>
inline real
FixedQuadrupleListTypesLambdaInteractionTemplate<_DihedralPotential>::
computeEnergy() {
  LOG4ESPP_INFO(theLogger, "compute energy of the quadruples");

  const bc::BC &bc = *getSystemRef().bc;  // boundary conditions
  real e = 0.0;
  for (FixedQuadrupleListLambda::IteratorParticleLambda it(fixedquadrupleList->getParticleQuadruples()); it.isValid(); ++it) {
    const Particle &p1 = *it->first;
    const Particle &p2 = *it->second;
    const Particle &p3 = *it->third;
    const Particle &p4 = *it->fourth;
    real lambda = it->lambda;

    longint type1 = p1.type();
    longint type2 = p2.type();
    longint type3 = p3.type();
    longint type4 = p4.type();

    const Potential &potential = getPotential(type1, type2, type3, type4);

    Real3D dist21, dist32, dist43;

    bc.getMinimumImageVectorBox(dist21, p2.position(), p1.position());
    bc.getMinimumImageVectorBox(dist32, p3.position(), p2.position());
    bc.getMinimumImageVectorBox(dist43, p4.position(), p3.position());

    e += lambda*potential.computeEnergy(dist21, dist32, dist43);
  }
  real esum;
  boost::mpi::all_reduce(*mpiWorld, e, esum, std::plus<real>());
  return esum;
}

template < typename _DihedralPotential > inline real
FixedQuadrupleListTypesLambdaInteractionTemplate < _DihedralPotential >::
computeEnergyDeriv() {
  std::cout << "Warning! At the moment computeEnergyDeriv() in _FixedQuadrupleListTypesLambdaInteractionTemplate does not work." << std::endl;
  return 0.0;
}

template<typename _DihedralPotential>
inline real
FixedQuadrupleListTypesLambdaInteractionTemplate<_DihedralPotential>::
computeEnergyAA() {
  return 0.0;
}

template<typename _DihedralPotential>
inline real
FixedQuadrupleListTypesLambdaInteractionTemplate<_DihedralPotential>::
computeEnergyCG() {
  return 0.0;
}

template<typename _DihedralPotential>
inline void
FixedQuadrupleListTypesLambdaInteractionTemplate<_DihedralPotential>::
computeVirialX(std::vector<real> &p_xx_total, int bins) {
}

template<typename _DihedralPotential>
inline real
FixedQuadrupleListTypesLambdaInteractionTemplate<_DihedralPotential>::
computeVirial() {
  LOG4ESPP_INFO(theLogger, "compute scalar virial of the quadruples");

  real w = 0.0;
  const bc::BC &bc = *getSystemRef().bc;  // boundary conditions
  for (FixedQuadrupleListLambda::IteratorParticleLambda it(fixedquadrupleList->getParticleQuadruples()); it.isValid(); ++it) {
    const Particle &p1 = *it->first;
    const Particle &p2 = *it->second;
    const Particle &p3 = *it->third;
    const Particle &p4 = *it->fourth;
    real lambda = it->lambda;

    longint type1 = p1.type();
    longint type2 = p2.type();
    longint type3 = p3.type();
    longint type4 = p4.type();

    const Potential &potential = getPotential(type1, type2, type3, type4);

    Real3D dist21, dist32, dist43;

    bc.getMinimumImageVectorBox(dist21, p2.position(), p1.position());
    bc.getMinimumImageVectorBox(dist32, p3.position(), p2.position());
    bc.getMinimumImageVectorBox(dist43, p4.position(), p3.position());

    Real3D force1, force2, force3, force4;

    potential.computeForce(force1, force2, force3, force4,
                           dist21, dist32, dist43);

    // TODO: formulas are not correct yet?

    w += dist21 * lambda*force1 + dist32 * lambda*force2;
  }

  real wsum;
  boost::mpi::all_reduce(*mpiWorld, w, wsum, std::plus<real>());
  return w;
}

template<typename _DihedralPotential>
inline void
FixedQuadrupleListTypesLambdaInteractionTemplate<_DihedralPotential>::
computeVirialTensor(Tensor &w) {
  LOG4ESPP_INFO(theLogger, "compute the virial tensor of the quadruples");

  Tensor wlocal(0.0);
  const bc::BC &bc = *getSystemRef().bc;

  for (FixedQuadrupleListLambda::IteratorParticleLambda it(fixedquadrupleList->getParticleQuadruples()); it.isValid(); ++it) {
    const Particle &p1 = *it->first;
    const Particle &p2 = *it->second;
    const Particle &p3 = *it->third;
    const Particle &p4 = *it->fourth;
    real lambda = it->lambda;

    longint type1 = p1.type();
    longint type2 = p2.type();
    longint type3 = p3.type();
    longint type4 = p4.type();

    const Potential &potential = getPotential(type1, type2, type3, type4);

    Real3D dist21, dist32, dist43;

    bc.getMinimumImageVectorBox(dist21, p2.position(), p1.position());
    bc.getMinimumImageVectorBox(dist32, p3.position(), p2.position());
    bc.getMinimumImageVectorBox(dist43, p4.position(), p3.position());

    Real3D force1, force2, force3, force4;

    potential.computeForce(force1, force2, force3, force4,
                           dist21, dist32, dist43);

    // TODO: formulas are not correct yet

    wlocal += Tensor(dist21, lambda*force1) - Tensor(dist32, lambda*force2);
  }
  // reduce over all CPUs
  Tensor wsum(0.0);
  boost::mpi::all_reduce(*mpiWorld, (double *) &wlocal, 6, (double *) &wsum, std::plus<double>());
  w += wsum;
}


template<typename _DihedralPotential>
inline void
FixedQuadrupleListTypesLambdaInteractionTemplate<_DihedralPotential>::
computeVirialTensor(Tensor &w, real z) {
  LOG4ESPP_INFO(theLogger, "compute the virial tensor of the quadruples");

  Tensor wlocal(0.0);
  const bc::BC &bc = *getSystemRef().bc;

  std::cout << "Warning!!! computeVirialTensor in specified volume doesn't work for "
      "_FixedQuadrupleListTypesLambdaInteractionTemplate at the moment" << std::endl;

  for (FixedQuadrupleListLambda::IteratorParticleLambda it(fixedquadrupleList->getParticleQuadruples()); it.isValid(); ++it) {
    const Particle &p1 = *it->first;
    const Particle &p2 = *it->second;
    const Particle &p3 = *it->third;
    const Particle &p4 = *it->fourth;
    real lambda = it->lambda;

    longint type1 = p1.type();
    longint type2 = p2.type();
    longint type3 = p3.type();
    longint type4 = p4.type();

    const Potential &potential = getPotential(type1, type2, type3, type4);

    Real3D dist21, dist32, dist43;

    bc.getMinimumImageVectorBox(dist21, p2.position(), p1.position());
    bc.getMinimumImageVectorBox(dist32, p3.position(), p2.position());
    bc.getMinimumImageVectorBox(dist43, p4.position(), p3.position());

    Real3D force1, force2, force3, force4;

    potential.computeForce(force1, force2, force3, force4,
                           dist21, dist32, dist43);

    // TODO: formulas are not correct yet

    wlocal += Tensor(dist21, lambda*force1) - Tensor(dist32, lambda*force2);
  }
  // reduce over all CPUs
  Tensor wsum(0.0);
  boost::mpi::all_reduce(*mpiWorld, (double *) &wlocal, 6, (double *) &wsum, std::plus<double>());
  w += wsum;
}

template<typename _DihedralPotential>
inline void
FixedQuadrupleListTypesLambdaInteractionTemplate<_DihedralPotential>::
computeVirialTensor(Tensor *w, int n) {
  std::cout << "Warning!!! computeVirialTensor in specified volume doesn't work for "
      "_FixedQuadrupleListTypesLambdaInteractionTemplate at the moment" << std::endl;
}

template<typename _DihedralPotential>
inline real
FixedQuadrupleListTypesLambdaInteractionTemplate<_DihedralPotential>::
getMaxCutoff() {
  real cutoff = 0.0;
  for (int i = 0; i < ntypes; i++) {
    for (int j = 0; j < ntypes; j++) {
      for (int k = 0; k < ntypes; k++) {
        for (int l = 0; l < ntypes; l++) {
          cutoff = std::max(cutoff, getPotential(i, j, k, l).getCutoff());
        }
      }
    }
  }
  return cutoff;
}
}  // namespace interaction
}  // namespace espressopp
#endif
