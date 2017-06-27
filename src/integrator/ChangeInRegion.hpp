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
#ifndef _INTEGRATOR_ChangeInRegion_HPP
#define _INTEGRATOR_ChangeInRegion_HPP

#include "types.hpp"
#include "logging.hpp"
#include "Extension.hpp"
#include "ParticleRegion.hpp"
#include "boost/signals2.hpp"
#include "Real3D.hpp"
#include "Particle.hpp"
#include "TopologyManager.hpp"

namespace espressopp {
namespace integrator {

class ChangeInRegion: public Extension {
 public:
  /*** Select particles with probability p. */
  ChangeInRegion(shared_ptr<System> system, shared_ptr<ParticleRegion> particle_region, real p);
  /*** Select only num_part */
  ChangeInRegion(shared_ptr<System> system, shared_ptr<ParticleRegion> particle_region, long num_part);
  /*** Select only percentage of particles that are in the region. */
  ChangeInRegion(shared_ptr<System> system, shared_ptr<ParticleRegion> particle_region, long xxx, real percentage);

  virtual ~ChangeInRegion() {
    disconnect();
  };

  void setParticleProperties(longint type_id, shared_ptr<TopologyParticleProperties> pp) {
    type_particleProperties.insert(std::make_pair(type_id, pp));
  }

  void setFlags(longint type_id, bool reset_velocity, bool reset_force, bool remove_particle) {

    int flag = 0;
    if (reset_velocity)
      flag |= R_VELOCITY;
    if (reset_force)
      flag |= R_FORCE;
    if (remove_particle)
      flag |= R_PARTICLE;

    LOG4ESPP_DEBUG(theLogger,
      "set flags type_id: " << type_id << " v:" << reset_velocity << " f:" << reset_force
      << " remove: " << remove_particle
      << " flag:" << flag);
    type_flags.insert(std::make_pair(type_id, flag));
  }

  /** Register this class so it can be used from Python. */
  static void registerPython();

 private:
  enum FLAGS {R_VELOCITY = 1, R_FORCE = 2, R_PARTICLE = 4};
  boost::signals2::connection sig_aftIntV;
  shared_ptr<ParticleRegion> particleRegion;
  real p_;  ///!< Probability of changing the values.
  longint num_particles_;
  real percentage_;
  std::string stats_filename_;

  std::map<longint, shared_ptr<TopologyParticleProperties> > type_particleProperties;
  std::map<longint, int> type_flags;

  void updateParticlesProb();
  void updateParticlesFlux();
  void updateParticlesPercentage();
  void connect();
  void disconnect();

  /** Logger */
  static LOG4ESPP_DECL_LOGGER(theLogger);
};

}
}

#endif
