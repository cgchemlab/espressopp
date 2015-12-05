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

class TopologyManager: public Extension {
 public:
  TopologyManager(shared_ptr<System> system);
  ~TopologyManager();

  /**
   * Pair list registration. The following types
   * are used to update those list based on the new
   * bond that appear in the system.
   */
  /**
   * Observe fixed pair list against changes.
   * @param fpl FixedPairList
   */
  void observeTuple(shared_ptr<FixedPairList> fpl,
                    longint type1, longint type2);
  /**
   * Register triple pair list with desire atom types.
   * @param ftl FixedTripleList
   * @param type1 The longint type_1
   * @param type2 The longint type_2
   * @param type3 The lontint type_3
   */
  void observeTriple(shared_ptr<FixedTripleList> ftl,
                     longint type1, longint type2, longint type3);
  /**
   * Register quadruple list with desire atom types.
   * @param fql FixedQuadrupleList
   * @param type1 The longint type_1
   * @param type2 The longint type_2
   * @param type3 The longint type_3
   * @param type4 The longint type_4
   */
  void observeQuadruple(shared_ptr<FixedQuadrupleList> fql,
                        longint type1, longint type2, longint type3, longint type4);

  /**
   * Initialized topology by looking for bonds in registered PairLists and
   * build adjacent list. Then this list is distributed among cpus so
   * everyone has the same adjacent list.
   */
  void InitializeTopology();
  /**
   * Print adjacent list.
   */
  void PrintTopology();

  /**
   * Get neighbour list.
   */
  python::list getNeighbourLists();

  /**
   * Rebuild the map of particle_id -> res_id and sync among cpus.
   */
  void Rebuild();

  static void registerPython();

  /// Defines map that stores res_id -> particle_id
  typedef std::set<longint> PSet;
  typedef std::map<longint, PSet *> ResParticleIds;

 private:
  typedef std::pair<longint, std::pair<longint, longint> > Triplets;
  typedef std::pair<longint, std::pair<longint, std::pair<longint, longint> > > Quadruplets;
  typedef std::vector<std::pair<longint, longint> > EdgesVector;

  /**
   * Handle signal from FixedPairList that new bond was created.
   */
  void onTupleAdded(longint pid1, longint pid2);
  /**
   * Process new bond.
   */
  void newBond(longint pid1, longint pid2);
  /**
   * Process new edge in topology.
   */
  void newEdge(longint pid1, longint pid2);
  /**
   * Update registered FixedTripleList with new entries.
   */
  void updateAngles(std::set<Triplets> &triplets);
  /**
   * Update registered FixedQuadrupleLists with new entries.
   */
  void updateDihedrals(std::set<Quadruplets> &quadruplets);
  /**
   * Based on topology, generate missing angles and dihedrals around newly
   * created bond.
   */
  void generateAnglesDihedrals(longint pid1,
                               longint pid2,
                               std::set<Quadruplets> &quadruplets,
                               std::set<Triplets> &triplets);
  /**
   * Exchange new topology and res_id data among cpus.
   */
  void exchangeData();

  /**
   * Merge two sets that holds mapping particle_id -> res_id
   * and update res_id of corresponding particles.
   */
  void mergeResIdSets(longint res_id_a, longint res_id_b);
  /**
   * Connecting/Disconnecting to signals.
   */
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
  typedef std::map<longint, std::set<int> *> GraphMap;
  bool update_angles_dihedrals;

  std::vector<shared_ptr<FixedPairList> > tupleMap_;
  TripleMap tripleMap_;
  QuadrupleMap quadrupleMap_;

  EdgesVector newEdges_;

  /** Adjacent list. */
  GraphMap *graph_;

  /** Logger */
  static LOG4ESPP_DECL_LOGGER(theLogger);
};

}  // end namespace integrator
}  // end namespace espressoppp
#endif //_INTEGRATOR_TOPOLOGYMANAGER_H