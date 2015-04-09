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
#ifndef _INTERACTION_VERLETLISTNONRECIPROCALINTERACTIONTEMPLATE_HPP
#define _INTERACTION_VERLETLISTNONRECIPROCALINTERACTIONTEMPLATE_HPP

#include <algorithm>
#include <vector>
#include "types.hpp"
#include "Interaction.hpp"
#include "Real3D.hpp"
#include "Tensor.hpp"
#include "Particle.hpp"
#include "VerletList.hpp"
#include "esutil/Array2D.hpp"
#include "bc/BC.hpp"

#include "storage/Storage.hpp"

/**
 * This class implements special type of non-bonded potential. The force is only exert on selected
 * particle type (defined by active_type) but the reciprocal force is not applied.
 *
 * The aim of this potential is to bring defined particles into local equilibrium but on the other
 * hand those particles should not effect other particles.
 *
 * Warning: This potential violate third Newton law and can produce unwanted results.
 */

namespace espressopp {
namespace interaction {
template < typename _Potential >
class VerletListNonReciprocalInteractionTemplate: public Interaction {
 protected:
  typedef _Potential Potential;

 public:
  VerletListNonReciprocalInteractionTemplate(shared_ptr<VerletList> _verletList, int _active_type)
      : verletList(_verletList), active_type(_active_type) {
    potentialArray = esutil::Array2D<Potential, esutil::enlarge>(0, 0, Potential());
    ntypes = 0;
    LOG4ESPP_WARN(_Potential::theLogger,
        "Warning! You are using non-reciprocal force with this potential"
        << " that will only acts on particles with type " << _active_type);
  }

  virtual ~VerletListNonReciprocalInteractionTemplate() {}

  void setVerletList(shared_ptr < VerletList > _verletList) {
    verletList = _verletList;
  }

  shared_ptr<VerletList> getVerletList() {
    return verletList;
  }

  void setPotential(int type1, int type2, const Potential &potential) {
    ntypes = std::max(ntypes, std::max(type1+1, type2+1));
    potentialArray.at(type1, type2) = potential;
    LOG4ESPP_INFO(_Potential::theLogger,
        "added potential for type1=" << type1 << " type2=" << type2);
    if (type1 != type2) {  // add potential in the other direction
       potentialArray.at(type2, type1) = potential;
       LOG4ESPP_INFO(_Potential::theLogger,
           "automatically added the same potential for type1=" << type2 << " type2=" << type1);
    }
  }

  // this is used in the innermost force-loop
  Potential &getPotential(int type1, int type2) {
    return potentialArray.at(type1, type2);
  }

  shared_ptr<Potential> getPotentialPtr(int type1, int type2) {
    return  make_shared<Potential>(potentialArray.at(type1, type2));
  }


  virtual void addForces();
  virtual real computeEnergy() { return 0.0;  }
  virtual real computeEnergyAA() { return 0.0; }
  virtual real computeEnergyCG() { return 0.0; }
  virtual void computeVirialX(std::vector<real> &p_xx_total, int bins);
  virtual real computeVirial() { return 0.0; }
  virtual void computeVirialTensor(Tensor& w);
  virtual void computeVirialTensor(Tensor& w, real z);
  virtual void computeVirialTensor(Tensor *w, int n);
  virtual real getMaxCutoff();
  virtual int bondType() { return Nonbonded; }

 protected:
  int ntypes;
  shared_ptr<VerletList> verletList;
  esutil::Array2D<Potential, esutil::enlarge> potentialArray;
  int active_type;  /// Type of particle that only will fell the force.
};


/** Inline implementation of templated methods. */
template < typename _Potential > inline void
VerletListNonReciprocalInteractionTemplate < _Potential >::addForces() {
  LOG4ESPP_DEBUG(_Potential::theLogger, "loop over verlet list pairs and add forces");
  
  for (PairList::Iterator it(verletList->getPairs()); it.isValid(); ++it) {
    Particle &p1 = *it->first;
    Particle &p2 = *it->second;
    int type1 = p1.type();
    int type2 = p2.type();
    const Potential &potential = getPotential(type1, type2);
    // shared_ptr<Potential> potential = getPotential(type1, type2);

    Real3D force(0.0);
    if (potential._computeForce(force, p1, p2)) {
      bool stat = false;
      if (type1 == active_type) {
        p1.force() += force;
        stat = true;
      } else if (type2 == active_type) {
        p2.force() -= force;
        stat = true;
      }
      if (stat)
        LOG4ESPP_TRACE(_Potential::theLogger, "id1=" << p1.id()
                       << " id2=" << p2.id() << " force=" << force);
    }
  }
}

template <typename _Potential>
inline void VerletListNonReciprocalInteractionTemplate < _Potential >::computeVirialX(
    std::vector<real> &p_xx_total, int bins) {
}

template <typename _Potential>
inline void VerletListNonReciprocalInteractionTemplate < _Potential >::computeVirialTensor(Tensor& w) {
}

template <typename _Potential>
inline void VerletListNonReciprocalInteractionTemplate < _Potential >::computeVirialTensor(Tensor& w, real z) {
}

template <typename _Potential>
inline void VerletListNonReciprocalInteractionTemplate < _Potential >::computeVirialTensor(Tensor *w, int n) {
}

template <typename _Potential>
inline real VerletListNonReciprocalInteractionTemplate< _Potential >::getMaxCutoff() {
  real cutoff = 0.0;
  for (int i = 0; i < ntypes; i++) {
    for (int j = 0; j < ntypes; j++) {
        cutoff = std::max(cutoff, getPotential(i, j).getCutoff());
    }
  }
  return cutoff;
}

}  // namespace interaction
}  // namespace espressopp
#endif
