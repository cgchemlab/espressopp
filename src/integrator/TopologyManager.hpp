/*
  Copyright (C) 2015
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

#ifndef _INTEGRATOR_TOPOLOGYMANAGER_H
#define _INTEGRATOR_TOPOLOGYMANAGER_H

#include <vector>

#include "types.hpp"
#include "logging.hpp"
#include "Extension.hpp"
#include "Particle.hpp"
#include "FixedPairList.hpp"
#include "System.hpp"

#include "boost/unordered_map.hpp"

namespace espressopp {
namespace integrator {

class TopologyManager : public Extension {
  /// Defines map that stores res_id -> particle_id
  typedef boost::unordered_multimap<longint, longint> ResParticleIds;
 public:
  TopologyManager(shared_ptr<System> system);

  void observeTuple(shared_ptr<FixedPairList> fpl);

  void rebuild();
 private:
  void onTupleAdded(longint pid1, longint pid2);
  void connect();
  void disconnect();

  shared_ptr<System> system_;
  ResParticleIds res_particle_ids_;
  bool is_dirty_;
  /** Logger */
  static LOG4ESPP_DECL_LOGGER(theLogger);
};

}  // end namespace integrator
}  // end namespace espressoppp
#endif //_INTEGRATOR_TOPOLOGYMANAGER_H
