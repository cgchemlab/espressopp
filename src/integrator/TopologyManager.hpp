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
  void PrintTopology();

  /** Rebuild map res_id -> particle_id */
  void Rebuild();

  static void registerPython();

  /// Defines map that stores res_id -> particle_id
  typedef std::set<longint> PSet;
  typedef std::map<longint, PSet*> ResParticleIds;

 private:
  typedef std::pair<longint, std::pair<longint, longint> > Triplets;
  typedef std::pair<longint, std::pair<longint, std::pair<longint, longint> > > Quadruplets;
  typedef std::vector<std::pair<longint, longint> > EdgesVector;

  /** Handle local tuple update. */
  void onTupleAdded(longint pid1, longint pid2);
  void newBond(longint pid1, longint pid2);
  void newEdge(longint pid1, longint pid2);
  void updateAngles(std::set<Triplets> &triplets);
  void updateDihedrals(std::set<Quadruplets> &quadruplets);
  void generateAnglesDihedrals(longint pid1, longint pid2, std::set<Quadruplets> &quadruplets, std::set<Triplets> &triplets);
  void exchangeData();
  void mergeResIdSets(longint res_id_a, longint res_id_b);
  void connect();
  void disconnect();

  shared_ptr<System> system_;
  ResParticleIds res_particle_ids_;
  std::vector<std::pair<longint, longint> > merge_sets_;

  boost::signals2::connection aftIntV2_, aftCalcF_;

  // Maping for tuples, triplets and quadruplets.
  typedef boost::unordered_map<
      longint,
      boost::unordered_map<
          longint,
          boost::unordered_map<
              longint,
              shared_ptr<FixedTripleList> > > > TripleMap;
  typedef boost::unordered_map<
      longint,
      boost::unordered_map<
          longint,
          boost::unordered_map<
              longint,
              boost::unordered_map<
                  longint,
                  shared_ptr<FixedQuadrupleList> > > > > QuadrupleMap;
  typedef std::map<longint, std::set<int>* > GraphMap;

  std::vector<shared_ptr<FixedPairList> > tupleMap_;
  TripleMap tripleMap_;
  QuadrupleMap quadrupleMap_;

  EdgesVector newEdges_;

  GraphMap *graph_;

  /** Logger */
  static LOG4ESPP_DECL_LOGGER(theLogger);
};

}  // end namespace integrator
}  // end namespace espressoppp
#endif //_INTEGRATOR_TOPOLOGYMANAGER_H
