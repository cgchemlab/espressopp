/*
  Copyright (c) 2015
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
#ifndef _INTERACTION_VERLETLISTINTERACTIONDYNAMICRESOLUTIONTEMPLATE_HPP
#define _INTERACTION_VERLETLISTINTERACTIONDYNAMICRESOLUTIONTEMPLATE_HPP

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
class VerletListDynamicResolutionInteractionTemplate: public Interaction {

protected:
  typedef _Potential Potential;

public:
  VerletListDynamicResolutionInteractionTemplate(
      shared_ptr<VerletList> _verletList, bool _cg_potential):
          verletList(_verletList), cgPotential(_cg_potential) {
    potentialArray = esutil::Array2D<Potential, esutil::enlarge>(0, 0, Potential());
    ntypes = 0;
    has_max_force_ = false;
  }

  virtual ~VerletListDynamicResolutionInteractionTemplate() {};

  void setVerletList(shared_ptr<VerletList> _verletList) {
    verletList = _verletList;
  }

  shared_ptr<VerletList> getVerletList() {
    return verletList;
  }

  void setPotential(int type1, int type2, const Potential &potential) {
    // typeX+1 because i<ntypes
    ntypes = std::max(ntypes, std::max(type1+1, type2+1));
    potentialArray.at(type1, type2) = potential;
    LOG4ESPP_INFO(
        _Potential::theLogger,
        "added potential for type1=" << type1 << " type2=" << type2);
    if (type1 != type2) { // add potential in the other direction
      potentialArray.at(type2, type1) = potential;
      LOG4ESPP_INFO(
          _Potential::theLogger,
          "automatically added the same potential for type1=" << type2 << " type2=" << type1);
    }
  }

  // this is used in the innermost force-loop
  Potential &getPotential(int type1, int type2) {
    return potentialArray.at(type1, type2);
  }

  // this is mainly used to access the potential from Python (e.g. to change parameters of the potential)
  shared_ptr<Potential> getPotentialPtr(int type1, int type2) {
    return  make_shared<Potential>(potentialArray.at(type1, type2));
  }

  void setMaxForce(real max_force) {
    max_force_ = max_force;
    has_max_force_ = true;
  }

  virtual void addForces();
  virtual real computeEnergy();
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
  bool cgPotential;
  real max_force_;
  bool has_max_force_;
};

//////////////////////////////////////////////////
// INLINE IMPLEMENTATION
//////////////////////////////////////////////////
template <typename _Potential> inline void
VerletListDynamicResolutionInteractionTemplate<_Potential>::
addForces() {
  LOG4ESPP_DEBUG(_Potential::theLogger, "loop over verlet list pairs and add forces");
  const bc::BC& bc = *(verletList->getSystemRef()).bc;  // boundary conditions
  for (PairList::Iterator it(verletList->getPairs()); it.isValid(); ++it) {
    Particle &p1 = *it->first;
    Particle &p2 = *it->second;
    int type1 = p1.type();
    int type2 = p2.type();
    real w12 = p1.lambda() * p2.lambda();
    if (cgPotential) {
      w12 = (1.0-w12);
    }
    if (w12 > 0.0) {
      const Potential &potential = getPotential(type1, type2);
      // shared_ptr<Potential> potential = getPotential(type1, type2);

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
        p1.force() += w12*force;
        p2.force() -= w12*force;
      }
    }
  }
}

template <typename _Potential>
inline real
VerletListDynamicResolutionInteractionTemplate<_Potential>::
computeEnergy() {
  LOG4ESPP_DEBUG(_Potential::theLogger, "loop over verlet list pairs and sum up potential energies");

  real es = 0.0;
  for (PairList::Iterator it(verletList->getPairs()); it.isValid(); ++it) {
    Particle &p1 = *it->first;
    Particle &p2 = *it->second;
    int type1 = p1.type();
    int type2 = p2.type();
    real w12 = p1.lambda() * p2.lambda();
    if (cgPotential) {
      LOG4ESPP_DEBUG(_Potential::theLogger, "cgPotential=" << cgPotential << " w12=" << w12);
      w12 = (1.0 - w12);
    }
    const Potential &potential = getPotential(type1, type2);

    real e = potential._computeEnergy(p1, p2);
    if (isnan(e) || isinf(e)) {
      e = 0.0;
    }
    es += w12*e;
  }

  // reduce over all CPUs
  real esum;
  boost::mpi::all_reduce(*getVerletList()->getSystem()->comm, es, esum, std::plus<real>());
  return esum;
}

template <typename _Potential> inline real
VerletListDynamicResolutionInteractionTemplate<_Potential>::
computeEnergyAA() {
  if (!cgPotential)
    return computeEnergy();
  return 0.0;
}

template < typename _Potential > inline real
VerletListDynamicResolutionInteractionTemplate < _Potential >::
computeEnergyCG() {
  if (cgPotential)
    return computeEnergy();
  return 0.0;
}

template <typename _Potential>
inline void
VerletListDynamicResolutionInteractionTemplate<_Potential>::
computeVirialX(std::vector<real> &p_xx_total, int bins) {
  LOG4ESPP_WARN(_Potential::theLogger, "Warning! computeVirialX() is not yet implemented.");
}

template <typename _Potential> inline real
VerletListDynamicResolutionInteractionTemplate<_Potential>::
computeVirial() {
  LOG4ESPP_DEBUG(_Potential::theLogger, "loop over verlet list pairs and sum up virial");
  const bc::BC& bc = *(verletList->getSystemRef()).bc;  // boundary conditions

  real w = 0.0;
  for (PairList::Iterator it(verletList->getPairs());
       it.isValid(); ++it) {
    Particle &p1 = *it->first;
    Particle &p2 = *it->second;
    int type1 = p1.type();
    int type2 = p2.type();
    real w12 = p1.lambda() * p2.lambda();
    if (cgPotential) {
      w12 = (1.0 - w12);
    }
    if (w12 > 0.0) {
      const Potential &potential = getPotential(type1, type2);
      // shared_ptr<Potential> potential = getPotential(type1, type2);

      Real3D force(0.0, 0.0, 0.0);
      Real3D dist;
      bc.getMinimumImageVectorBox(dist, p1.position(), p2.position());
      if (potential._computeForce(force, p1, p2)) {
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
        Real3D r21 = p1.position() - p2.position();
        w = w + w12 * r21 * force;
      }
    }
  }

  // reduce over all CPUs
  real wsum;
  boost::mpi::all_reduce(*mpiWorld, w, wsum, std::plus<real>());
  return wsum;
}

template <typename _Potential> inline void
VerletListDynamicResolutionInteractionTemplate<_Potential>::
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
    real w12 = p1.lambda() * p2.lambda();
    if (cgPotential) {
      w12 = (1.0 - w12);
    }
    if (w12 > 0.0) {
      const Potential &potential = getPotential(type1, type2);
      // shared_ptr<Potential> potential = getPotential(type1, type2);

      Real3D force(0.0, 0.0, 0.0);
      Real3D dist;
      bc.getMinimumImageVectorBox(dist, p1.position(), p2.position());
      if (potential._computeForce(force, p1, p2)) {
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
        wlocal += Tensor(dist, w12 * force);
      }
    }
  }

  // reduce over all CPUs
  Tensor wsum(0.0);
  boost::mpi::all_reduce(*mpiWorld, (double*)&wlocal, 6, (double*)&wsum, std::plus<double>());
  w += wsum;
}

// local pressure tensor for layer, plane is defined by z coordinate
template <typename _Potential> inline void
VerletListDynamicResolutionInteractionTemplate<_Potential>::
computeVirialTensor(Tensor& w, real z) {
  LOG4ESPP_ERROR(_Potential::theLogger, "Not implemented!");
}

// it will calculate the pressure in 'n' layers along Z axis
// the first layer has coordinate 0.0 the last - (Lz - Lz/n)
template <typename _Potential> inline void
VerletListDynamicResolutionInteractionTemplate<_Potential>::
computeVirialTensor(Tensor *w, int n) {
  LOG4ESPP_ERROR(_Potential::theLogger, "Not implemented!");
}

template <typename _Potential>
inline real VerletListDynamicResolutionInteractionTemplate<_Potential>::
getMaxCutoff() {
  real cutoff = 0.0;
  for (int i = 0; i < ntypes; i++) {
    for (int j = 0; j < ntypes; j++) {
      cutoff = std::max(cutoff, getPotential(i, j).getCutoff());
      // cutoff = std::max(cutoff, getPotential(i, j)->getCutoff());
    }
  }
  return cutoff;
}
}  // end namespace interaction
}  // end namespace espressopp
#endif
