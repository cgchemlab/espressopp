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

#ifndef _INTEGRATOR_ATRPActivator_HPP
#define _INTEGRATOR_ATRPActivator_HPP

#include <algorithm>
#include <esutil/Timer.hpp>
#include "types.hpp"
#include "logging.hpp"
#include "Extension.hpp"
#include "boost/signals2.hpp"
#include "Real3D.hpp"
#include "Particle.hpp"
#include "iterator/CellListIterator.hpp"
#include "storage/NodeGrid.hpp"
#include "storage/DomainDecomposition.hpp"

namespace espressopp {
namespace integrator {

struct ReactiveCenter {
  longint state;  ///<! minimal chemical state
  bool is_activator;
  longint delta_state;  ///<! update of chemical potential
  shared_ptr<ParticleProperties> new_property;  ///<! new property

  ReactiveCenter() {
    state = -1;
    delta_state = 0;
  }

  ReactiveCenter(longint state_, bool is_activator_, longint delta_state_,
                 shared_ptr<ParticleProperties> new_property_) {
    state = state_;
    is_activator = is_activator_;
    delta_state = delta_state_;
    new_property = new_property_;
  }
};

class ATRPActivator: public Extension {
 public:
  ATRPActivator(shared_ptr<System> system, longint interval, longint num_per_interval, real ratio_activator,
                real ratio_deactivator, real delta_catalyst, real k_activate, real k_deactivate);

  virtual ~ATRPActivator() {};

  void addReactiveCenter(longint type_id,
                         longint state,
                         bool is_activator,
                         shared_ptr<ParticleProperties> pp,
                         longint delta_state);

  /** Register this class so it can be used from Python. */
  static void registerPython();

 private:
  boost::signals2::connection sig_aftIntV;
  longint interval_;
  longint num_per_interval_;
  real ratio_activator_;
  real ratio_deactivator_;
  real delta_catalyst_;
  real k_activate_;
  real k_deactivate_;
  std::string stats_filename_;

  std::vector<real> stats_k_activator;

  typedef boost::unordered_multimap<longint, ReactiveCenter> SpeciesMap;
  SpeciesMap species_map_;

  shared_ptr<esutil::RNG> rng_;

  void updateParticles();
  void updateGhost(const std::vector<Particle *> &modified_particles);

  void saveStatistics(std::string filename);

  void connect();
  void disconnect();

  /** Logger */
  static LOG4ESPP_DECL_LOGGER(theLogger);

  /** Timers */
  real timeUpdateParticles;
  real timeUpdateGhost;
  esutil::WallTimer wallTimer;

  void resetTimers() {
      wallTimer.reset();

      timeUpdateParticles = 0.0;
      timeUpdateGhost = 0.0;
  }

  python::list getTimers();
};

}
}

#endif
