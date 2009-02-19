#ifndef _INTERACTION_INTERACTION_HPP
#define _INTERACTION_INTERACTION_HPP

//base class

#include "types.hpp"
#include "particlestorage/ParticleStorage.hpp"
#include "particleset/ParticleSet.hpp"

namespace espresso {
  namespace interaction {
    class Interaction {
    public:
      virtual ~Interaction() {}
      virtual real computeEnergy(const Real3D &dist, 
				 const espresso::particleset::ParticleSet::const_reference p1,
				 const espresso::particleset::ParticleSet::const_reference p2) const = 0;
      virtual Real3D computeForce(const Real3D &dist, 
				  const espresso::particleset::ParticleSet::const_reference p1,
				  const espresso::particleset::ParticleSet::const_reference p2) const = 0;
      virtual real getCutoff() const = 0;
      virtual real getCutoffSqr() const = 0;
    };
  }
}

#endif