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
#include <map>

#include "types.hpp"
#include "logging.hpp"
#include "Extension.hpp"
#include "Particle.hpp"
#include "FixedPairList.hpp"
#include "System.hpp"

namespace espressopp {
namespace integrator {

class TopologyManager : public Extension {
 public:
  TopologyManager(shared_ptr<System> system);
  ~TopologyManager();

  /**
   * Observe fixed pair list against chages.
   * @param fpl FixedPairList
   */
  void observeTuple(shared_ptr<FixedPairList> fpl);

  /** Rebuild map res_id -> particle_id */
  void Rebuild();

  static void registerPython();

  /// Defines map that stores res_id -> particle_id
  typedef std::set<longint> PSet;
  typedef std::map<longint, PSet*> ResParticleIds;

 private:
  /** Handle local tuple update. */
  void onTupleAdded(longint pid1, longint pid2);
  void exchangeData();
  void mergeResIdSets(longint res_id_a, longint res_id_b);
  void connect();
  void disconnect();

  shared_ptr<System> system_;
  ResParticleIds res_particle_ids_;
  std::vector<longint> merge_sets_;

  boost::signals2::connection aftIntV_;

  /** Logger */
  static LOG4ESPP_DECL_LOGGER(theLogger);
};

}  // end namespace integrator
}  // end namespace espressoppp
#endif //_INTEGRATOR_TOPOLOGYMANAGER_H
