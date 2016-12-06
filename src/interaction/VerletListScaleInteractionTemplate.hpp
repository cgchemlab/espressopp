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
#ifndef _INTERACTION_VerletListScaleInteractionTemplate_HPP
#define _INTERACTION_VerletListScaleInteractionTemplate_HPP

#include <esutil/ParticlePairScaling.hpp>
#include "types.hpp"
#include "Interaction.hpp"
#include "Real3D.hpp"
#include "Tensor.hpp"
#include "Particle.hpp"
#include "VerletList.hpp"
#include "esutil/Array2D.hpp"
#include "bc/BC.hpp"

#include "storage/Storage.hpp"

namespace espressopp {
namespace interaction {

template <typename _Potential>
class VerletListScaleInteractionTemplate: public Interaction {

protected:
  typedef _Potential Potential;

public:
  VerletListScaleInteractionTemplate(shared_ptr<VerletList> _verletList, shared_ptr<esutil::ParticlePairScaling> pps)
      : verletList(_verletList), pair_scaling_(pps) {
    potentialArray = esutil::Array2D<Potential, esutil::enlarge>(0, 0, Potential());
    ntypes = 0;
    has_max_force_ = false;
  }

  virtual ~VerletListScaleInteractionTemplate() {};

  void setVerletList(shared_ptr<VerletList> _verletList) {
    verletList = _verletList;
  }

  shared_ptr<VerletList> getVerletList() {
    return verletList;
  }

  void setPotential(int type1, int type2, const Potential &potential) {
    ntypes = std::max(ntypes, std::max(type1+1, type2+1));
    potentialArray.at(type1, type2) = potential;
    LOG4ESPP_INFO(
        _Potential::theLogger, "added potential for type1=" << type1 << " type2=" << type2);
    if (type1 != type2) { // add potential in the other direction
      potentialArray.at(type2, type1) = potential;
      LOG4ESPP_INFO(
          _Potential::theLogger, "automatically added the same potential for type1=" << type2 << " type2=" << type1);
    }
  }

  Potential &getPotential(int type1, int type2) {
    return potentialArray.at(type1, type2);
  }

  shared_ptr<Potential> getPotentialPtr(int type1, int type2) {
    return  make_shared<Potential>(potentialArray.at(type1, type2));
  }

  void setMaxForce(real max_force) {
    max_force_ = max_force;
    has_max_force_ = true;
  }

  virtual void addForces();
  virtual real computeEnergy();
  virtual real computeEnergyDeriv();
  virtual real computeEnergyAA();
  virtual real computeEnergyCG();
  virtual void computeVirialX(std::vector<real> &p_xx_total, int bins);
  virtual real computeVirial();
  virtual void computeVirialTensor(Tensor& w);
  virtual void computeVirialTensor(Tensor& w, real z);
  virtual void computeVirialTensor(Tensor *w, int n);
  virtual real getMaxCutoff();
  virtual int bondType() { return Nonbonded; }

protected:
  int ntypes;
  shared_ptr<VerletList> verletList;
  esutil::Array2D<Potential, esutil::enlarge> potentialArray;
  real max_force_;
  bool has_max_force_;
  shared_ptr<esutil::ParticlePairScaling> pair_scaling_;
};

template <typename _Potential> inline void
VerletListScaleInteractionTemplate<_Potential>::
addForces() {
  LOG4ESPP_DEBUG(_Potential::theLogger, "loop over verlet list pairs and add forces");
  const bc::BC& bc = *(verletList->getSystemRef()).bc;  // boundary conditions
  for (PairList::Iterator it(verletList->getPairs()); it.isValid(); ++it) {
    Particle &p1 = *it->first;
    Particle &p2 = *it->second;
    int type1 = p1.type();
    int type2 = p2.type();
    const Potential &potential = getPotential(type1, type2);

    Real3D force(0.0);
    Real3D dist;
    bc.getMinimumImageVectorBox(dist, p1.position(), p2.position());
    if(potential._computeForce(force, p1, p2, dist)) {
      if (has_max_force_) {
        if (force.isNaNInf()) {
          force = (dist / dist.abs()) * max_force_;
        } else {
          real abs_force = force.abs();
          if (abs_force > max_force_) {
            force = (force / abs_force) * max_force_;
          }
        }
      }
      real pair_scaling = pair_scaling_->getPairScaling(p1.id(), p2.id());
      p1.force() += pair_scaling*force;
      p2.force() -= pair_scaling*force;
    }
  }
}

template <typename _Potential>
inline real
VerletListScaleInteractionTemplate<_Potential>::
computeEnergy() {
  LOG4ESPP_DEBUG(_Potential::theLogger, "loop over verlet list pairs and sum up potential energies");

  real es = 0.0;
  for (PairList::Iterator it(verletList->getPairs()); it.isValid(); ++it) {
    Particle &p1 = *it->first;
    Particle &p2 = *it->second;
    int type1 = p1.type();
    int type2 = p2.type();

    const Potential &potential = getPotential(type1, type2);

    real e = potential._computeEnergy(p1, p2);
    if (isnan(e) || isinf(e)) {
      e = 0.0;
    } else {
      real pair_scaling = pair_scaling_->getPairScaling(p1.id(), p2.id());
      e *= pair_scaling;
    }
    es += e;
  }

  // reduce over all CPUs
  real esum;
  boost::mpi::all_reduce(*getVerletList()->getSystem()->comm, es, esum, std::plus<real>());
  return esum;
}

template <typename _Potential> inline real
VerletListScaleInteractionTemplate<_Potential>::
computeEnergyAA() {
  return computeEnergy();
}

template < typename _Potential > inline real
VerletListScaleInteractionTemplate < _Potential >::
computeEnergyCG() {
  return 0.0;
}
template < typename _Potential > inline real
VerletListScaleInteractionTemplate < _Potential >::
computeEnergyDeriv() {
  LOG4ESPP_WARN(_Potential::theLogger, "Warning! computeEnergyDeriv() is not yet implemented.");
  return 0.0;
}

template <typename _Potential>
inline void
VerletListScaleInteractionTemplate<_Potential>::
computeVirialX(std::vector<real> &p_xx_total, int bins) {
  LOG4ESPP_WARN(_Potential::theLogger, "Warning! computeVirialX() is not yet implemented.");
}

template <typename _Potential> inline real
VerletListScaleInteractionTemplate<_Potential>::
computeVirial() {
  LOG4ESPP_DEBUG(_Potential::theLogger, "loop over verlet list pairs and sum up virial");
  const bc::BC& bc = *(verletList->getSystemRef()).bc;  // boundary conditions

  real w = 0.0;
  for (PairList::Iterator it(verletList->getPairs()); it.isValid(); ++it) {
    Particle &p1 = *it->first;
    Particle &p2 = *it->second;
    int type1 = p1.type();
    int type2 = p2.type();
    const Potential &potential = getPotential(type1, type2);

    Real3D force(0.0, 0.0, 0.0);
    Real3D dist;
    bc.getMinimumImageVectorBox(dist, p1.position(), p2.position());
    if (potential._computeForce(force, p1, p2, dist)) {
      if (has_max_force_) {
        if (force.isNaNInf()) {
          force = (dist/dist.abs()) * max_force_;
        } else {
          real abs_force = force.abs();
          if (abs_force > max_force_) {
            force = (force/abs_force) * max_force_;
          }
        }
      }
      real pair_scaling = pair_scaling_->getPairScaling(p1.id(), p2.id());
      Real3D r21 = p1.position() - p2.position();
      w = w + pair_scaling * r21 * force;
    }
  }

  real wsum;
  boost::mpi::all_reduce(*mpiWorld, w, wsum, std::plus<real>());
  return wsum;
}

template <typename _Potential> inline void
VerletListScaleInteractionTemplate<_Potential>::
computeVirialTensor(Tensor& w) {
  LOG4ESPP_DEBUG(_Potential::theLogger, "loop over verlet list pairs and sum up virial tensor");
  const bc::BC& bc = *(verletList->getSystemRef()).bc;  // boundary conditions
  Tensor wlocal(0.0);
  for (PairList::Iterator it(verletList->getPairs());
       it.isValid(); ++it) {
    Particle &p1 = *it->first;
    Particle &p2 = *it->second;
    int type1 = p1.type();
    int type2 = p2.type();
    const Potential &potential = getPotential(type1, type2);
    Real3D force(0.0, 0.0, 0.0);
    Real3D dist;
    bc.getMinimumImageVectorBox(dist, p1.position(), p2.position());
    if (potential._computeForce(force, p1, p2, dist)) {
      if (has_max_force_) {
        if (force.isNaNInf()) {
          force = (dist/dist.abs()) * max_force_;
        } else {
          real abs_force = force.abs();
          if (abs_force > max_force_) {
            force = (force/abs_force) * max_force_;
          }
        }
      }
      real pair_scaling = pair_scaling_->getPairScaling(p1.id(), p2.id());
      wlocal += Tensor(dist, pair_scaling*force);
    }
  }

  // reduce over all CPUs
  Tensor wsum(0.0);
  boost::mpi::all_reduce(*mpiWorld, (double*)&wlocal, 6, (double*)&wsum, std::plus<double>());
  w += wsum;
}

template <typename _Potential> inline void
VerletListScaleInteractionTemplate<_Potential>::
computeVirialTensor(Tensor& w, real z) {
  LOG4ESPP_ERROR(_Potential::theLogger, "Not implemented!");
}

template <typename _Potential> inline void
VerletListScaleInteractionTemplate<_Potential>::
computeVirialTensor(Tensor *w, int n) {
  LOG4ESPP_ERROR(_Potential::theLogger, "Not implemented!");
}

template <typename _Potential>
inline real VerletListScaleInteractionTemplate<_Potential>::
getMaxCutoff() {
  real cutoff = 0.0;
  for (int i = 0; i < ntypes; i++) {
    for (int j = 0; j < ntypes; j++) {
      cutoff = std::max(cutoff, getPotential(i, j).getCutoff());
    }
  }
  return cutoff;
}
}  // end namespace interaction
}  // end namespace espressopp
#endif
