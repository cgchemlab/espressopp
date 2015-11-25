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
#include "FixedTripleList.hpp"
#include "FixedQuadrupleList.hpp"
#include "System.hpp"

namespace espressopp {
namespace integrator {

class TopologyManager : public Extension {
 public:
  TopologyManager(shared_ptr<System> system);
  ~TopologyManager();

  /**
   * Observe fixed pair list against changes.
   * @param fpl FixedPairList
   */
  void observeTuple(shared_ptr<FixedPairList> fpl,
                    longint type1, longint type2);
  void observeTriple(shared_ptr<FixedTripleList> ftl,
                     longint type1, longint type2, longint type3);
  void observeQuadruple(shared_ptr<FixedQuadrupleList> fql,
                        longint type1, longint type2, longint type3, longint type4);

  /**
   * Initial topology.
   */
  void InitializeTopology();

  /** Rebuild map res_id -> particle_id */
  void Rebuild();

  static void registerPython();

  /// Defines map that stores res_id -> particle_id
  typedef std::set<longint> PSet;
  typedef std::map<longint, PSet*> ResParticleIds;

 private:
  /** Handle local tuple update. */
  void onTupleAdded(longint pid1, longint pid2);
  void newBond(Particle *p1, Particle *p2);
  void exchangeData();
  void mergeResIdSets(longint res_id_a, longint res_id_b);
  void connect();
  void disconnect();

  shared_ptr<System> system_;
  ResParticleIds res_particle_ids_;
  std::vector<longint> merge_sets_;

  boost::signals2::connection aftIntV_;

  // Maping for tuples, triplets and quadruplets.
  typedef boost::unordered_map<std::pair<longint, longint>, shared_ptr<FixedPairList> > TupleMap;
  typedef boost::unordered_map<boost::tuple<longint, longint, longint>,
                               shared_ptr<FixedTripleList> > TripleMap;
  typedef boost::unordered_map<boost::tuple<longint, longint, longint, longint>,
                               shared_ptr<FixedQuadrupleList> > QuadrupleMap;

  TupleMap tupleMap_;
  TripleMap tripleMap_;
  QuadrupleMap quadrupleMap_;

  /** Logger */
  static LOG4ESPP_DECL_LOGGER(theLogger);
};

}  // end namespace integrator
}  // end namespace espressoppp
#endif //_INTEGRATOR_TOPOLOGYMANAGER_H
