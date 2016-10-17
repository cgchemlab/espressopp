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
#ifndef _INTEGRATOR_ChangeParticleType_HPP
#define _INTEGRATOR_ChangeParticleType_HPP

#include "types.hpp"
#include "logging.hpp"
#include "Extension.hpp"
#include "ParticleGroup.hpp"
#include "boost/signals2.hpp"
#include "Real3D.hpp"
#include "Particle.hpp"
#include "storage/NodeGrid.hpp"
#include "storage/DomainDecomposition.hpp"

namespace espressopp {
namespace integrator {

class ChangeParticleType: public Extension {
 public:
  ChangeParticleType(
      shared_ptr<System> system,
      longint interval,
      longint num_per_interval,
      longint old_type,
      longint new_type);

  virtual ~ChangeParticleType() {};

  /** Register this class so it can be used from Python. */
  static void registerPython();

 private:
  boost::signals2::connection sig_aftIntV;
  shared_ptr<ParticleGroup> particle_group_;
  longint interval_;
  longint num_per_interval_;

  longint old_type_;
  longint new_type_;

  void updateParticles();
  void connect();
  void disconnect();

  /** Logger */
  static LOG4ESPP_DECL_LOGGER(theLogger);
  void updateGhost(const std::set<Particle *> &modified_particles);
};

}
}

#endif
