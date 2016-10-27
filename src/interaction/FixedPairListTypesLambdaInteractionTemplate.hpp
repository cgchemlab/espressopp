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

#ifndef _INTERACTION_FIXEDPAIRLAMBDALISTTYPESINTERACTIONTEMPLATE_HPP
#define _INTERACTION_FIXEDPAIRLAMBDALISTTYPESINTERACTIONTEMPLATE_HPP

#include "types.hpp"
#include "Interaction.hpp"
#include "Real3D.hpp"
#include "Tensor.hpp"
#include "Particle.hpp"
#include "FixedPairListLambda.hpp"
#include "esutil/Array2D.hpp"
#include "bc/BC.hpp"
#include "SystemAccess.hpp"
#include "types.hpp"

#include "storage/Storage.hpp"

namespace espressopp {
namespace interaction {
template<typename _Potential>
class FixedPairListTypesLambdaInteractionTemplate: public Interaction, SystemAccess {

 protected:
  typedef _Potential Potential;

 public:
  FixedPairListTypesLambdaInteractionTemplate
      (shared_ptr <System> system,
       shared_ptr <FixedPairListLambda> _fixedpairList)
      : SystemAccess(system), fixedpairList(_fixedpairList) {
    potentialArray = esutil::Array2D<Potential, esutil::enlarge>(0, 0, Potential());
    ntypes = 0;
  }

  virtual ~FixedPairListTypesLambdaInteractionTemplate() {};

  void
  setFixedPairList(shared_ptr <FixedPairListLambda> _fixedpairList) {
    fixedpairList = _fixedpairList;
  }

  shared_ptr <FixedPairListLambda> getFixedPairList() {
    return fixedpairList;
  }

  void
  setPotential(int type1, int type2, const Potential &potential) {
    ntypes = std::max(ntypes, std::max(type1 + 1, type2 + 1));
    potentialArray.at(type1, type2) = potential;
    if (type1 != type2) { // add potential in the other direction
      potentialArray.at(type2, type1) = potential;
    }
  }

  // this is used in the innermost force-loop
  Potential &getPotential(int type1, int type2) {
    return potentialArray.at(type1, type2);
  }

  // this is mainly used to access the potential from Python (e.g. to change parameters of the potential)
  shared_ptr <Potential> getPotentialPtr(int type1, int type2) {
    return make_shared<Potential>(potentialArray.at(type1, type2));
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
  esutil::Array2D <Potential, esutil::enlarge> potentialArray;
  esutil::Array2D <shared_ptr<Potential>, esutil::enlarge> potentialArrayPtr;
};

//////////////////////////////////////////////////
// INLINE IMPLEMENTATION
//////////////////////////////////////////////////
template<typename _Potential>
inline void
FixedPairListTypesLambdaInteractionTemplate<_Potential>::addForces() {
  LOG4ESPP_INFO(theLogger, "add forces computed by the FixedPair List");
  const bc::BC &bc = *getSystemRef().bc;

  for (FixedPairListLambda::Iterator it(*fixedpairList); it.isValid(); ++it) {
    Particle &p1 = *it->first;
    Particle &p2 = *it->second;
    real lambda = it->third;
    int type1 = p1.type();
    int type2 = p2.type();
    const Potential &potential = getPotential(type1, type2);

    Real3D force(0.0);

    Real3D dist;
    bc.getMinimumImageVectorBox(dist, p1.position(), p2.position());
    if (potential._computeForce(force, p1, p2, dist)) {
      p1.force() += lambda*force;
      p2.force() -= lambda*force;
    }
  }
}

template<typename _Potential>
inline real
FixedPairListTypesLambdaInteractionTemplate<_Potential>::
computeEnergy() {
  LOG4ESPP_INFO(theLogger, "compute energy of the FixedPair list pairs");

  real e = 0.0;
  real es = 0.0;
  const bc::BC &bc = *getSystemRef().bc;  // boundary conditions
  for (FixedPairListLambda::Iterator it(*fixedpairList); it.isValid(); ++it) {
    Particle &p1 = *it->first;
    Particle &p2 = *it->second;
    real lambda = it->third;
    int type1 = p1.type();
    int type2 = p2.type();
    const Potential &potential = getPotential(type1, type2);
    Real3D r21;
    bc.getMinimumImageVectorBox(r21, p1.position(), p2.position());
    e = lambda*potential._computeEnergy(p1, p2, r21);
    es += e;
    LOG4ESPP_TRACE(theLogger, "id1=" << p1.id() << " id2=" << p2.id() << " potential energy=" << e);
  }

  // reduce over all CPUs
  real esum;
  boost::mpi::all_reduce(*mpiWorld, es, esum, std::plus<real>());
  return esum;
}

template<typename _Potential>
inline real
FixedPairListTypesLambdaInteractionTemplate<_Potential>::
computeEnergyDeriv() {
  LOG4ESPP_ERROR(theLogger,
                 "Warning! At the moment computeEnergyDeriv() in "
                     << "FixedPairListTypesLambdaInteractionTemplate does not work.");
  return 0.0;
}

template<typename _Potential>
inline real
FixedPairListTypesLambdaInteractionTemplate<_Potential>::
computeEnergyAA() {
  LOG4ESPP_ERROR(theLogger,
                 "Warning! At the moment computeEnergyAA() in "
                     << "FixedPairListTypesLambdaInteractionTemplate does not work.");
  return 0.0;
}

template<typename _Potential>
inline real
FixedPairListTypesLambdaInteractionTemplate<_Potential>::
computeEnergyCG() {
  LOG4ESPP_ERROR(theLogger,
                 "Warning! At the moment computeEnergyCG() in "
                     << "FixedPairListTypesLambdaInteractionTemplate does not work.");
  return 0.0;
}

template<typename _Potential>
inline void
FixedPairListTypesLambdaInteractionTemplate<_Potential>::
computeVirialX(std::vector <real> &p_xx_total, int bins) {
  LOG4ESPP_ERROR(theLogger,
                 "Warning! At the moment computeVirialX() in "
                     << "FixedPairListTypesLambdaInteractionTemplate does not work.");
}

template<typename _Potential>
inline real
FixedPairListTypesLambdaInteractionTemplate<_Potential>::
computeVirial() {
  LOG4ESPP_INFO(theLogger, "compute the virial for the Fixed Pair List with types");

  real w = 0.0;
  const bc::BC &bc = *getSystemRef().bc;  // boundary conditions
  for (FixedPairListLambda::Iterator it(*fixedpairList);
       it.isValid(); ++it) {
    Particle &p1 = *it->first;
    Particle &p2 = *it->second;
    int type1 = p1.type();
    int type2 = p2.type();
    real lambda = it->third;
    const Potential &potential = getPotential(type1, type2);
    // shared_ptr<Potential> potential = getPotential(type1, type2);

    Real3D force(0.0, 0.0, 0.0);
    Real3D r21;
    bc.getMinimumImageVectorBox(r21, p1.position(), p2.position());
    if (potential._computeForce(force, p1, p2, r21)) {
      w = w + lambda*r21*force;
    }
  }

  // reduce over all CPUs
  real wsum;
  boost::mpi::all_reduce(*mpiWorld, w, wsum, std::plus<real>());
  return wsum;
}

template<typename _Potential>
inline void
FixedPairListTypesLambdaInteractionTemplate<_Potential>::computeVirialTensor(Tensor &w) {
  LOG4ESPP_ERROR(theLogger,
                 "Warning! At the moment computeVirialTensor() in "
                     << "FixedPairListTypesLambdaInteractionTemplate does not work.");
}

template<typename _Potential>
inline void
FixedPairListTypesLambdaInteractionTemplate<_Potential>::
computeVirialTensor(Tensor &w, real z) {
  LOG4ESPP_ERROR(theLogger,
                 "Warning! At the moment computeVirialTensor() in "
                     << "FixedPairListTypesLambdaInteractionTemplate does not work.");

}

template<typename _Potential>
inline void
FixedPairListTypesLambdaInteractionTemplate<_Potential>::
computeVirialTensor(Tensor *w, int n) {
  LOG4ESPP_ERROR(theLogger,
                 "Warning! At the moment computeVirialTensor() in "
                     << "FixedPairListTypesLambdaInteractionTemplate does not work.");
}

template<typename _Potential>
inline real
FixedPairListTypesLambdaInteractionTemplate<_Potential>::getMaxCutoff() {
  real cutoff = 0.0;
  for (int i = 0; i < ntypes; i++) {
    for (int j = 0; j < ntypes; j++) {
      cutoff = std::max(cutoff, getPotential(i, j).getCutoff());
    }
  }
  return cutoff;
}
}
}
#endif
