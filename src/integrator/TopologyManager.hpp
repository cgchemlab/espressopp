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
#include "FixedPairListLambda.hpp"
#include "FixedTripleList.hpp"
#include "FixedTripleListLambda.hpp"
#include "FixedQuadrupleList.hpp"
#include "FixedQuadrupleListLambda.hpp"
#include "System.hpp"
#include "esutil/Timer.hpp"
#include "boost/unordered_set.hpp"

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

  void observeTuple(shared_ptr<FixedPairList> fpl);

  /**
   * Observe fixed pair list against changes.
   * @param fpl FixedPairList
   * @param type1 particle type
   * @param type2 particle type
   * @param level tuple level; number of bonds between two particles on the list
   */
  void registerTuple(shared_ptr<FixedPairList> fpl, longint type1, longint type2);

  void register14Tuple(shared_ptr<FixedPairList> fpl, longint type1, longint type2);

  /**
   * Register triple pair list with desire atom types.
   * @param ftl FixedTripleList
   * @param type1 The longint type_1
   * @param type2 The longint type_2
   * @param type3 The lontint type_3
   */
  void registerTriple(shared_ptr<FixedTripleList> ftl,
                      longint type1, longint type2, longint type3);
  /**
   * Register quadruple list with desire atom types.
   * @param fql FixedQuadrupleList
   * @param type1 The longint type_1
   * @param type2 The longint type_2
   * @param type3 The longint type_3
   * @param type4 The longint type_4
   */
  void registerQuadruple(shared_ptr<FixedQuadrupleList> fql,
                         longint type1, longint type2, longint type3, longint type4);

  /**
   * Register the action to change neighbour property separated by nb_level from the root.
   *
   * @param type_id The type of particle that is neighbour of root separated by nb_level edges.
   * @param pp The new particle properties
   */
  void registerNeighbourPropertyChange(longint type_id, shared_ptr<ParticleProperties> pp, longint nb_level);

  /**
   * Register the action to remove bond that is `nb_level` bonds from the root particle
   * and involves particles of pid1, pid2.
   */
  void registerNeighbourBondToRemove(longint type_id, longint nb_level, longint type_pid1, longint type_pid2);

  /**
   * Interface for invoking property change of neighbour particles.
   *
   * @param root The root particle.
   */
  void invokeNeighbourPropertyChange(Particle &root);

  void invokeNeighbourBondRemove(Particle &root);

  /**
   * Schedule particle properties change. This is done after all other operations on topology. This is a local
   * operation, no communication is required.
   *
   * @param pid particle to change
   */
  void invokeParticlePropertiesChange(longint pid);
  void registerLocalPropertyChange(longint type_id, shared_ptr<ParticleProperties> pp);

  /**
   * Interface for checking if two residues are connected.
   * @param rid1 residue id
   * @param rid2 residue id
   * @return true if connected otherwise false
   */
  bool isResiduesConnected(longint rid1, longint rid2);

  /**
   * Check if two particles are connected.
   * @param pid1 particle id
   * @param pid2 particle id
   * @return true if connected otherwise false
   */
  bool isParticleConnected(longint pid1, longint pid2);

  bool isNeighbourParticleInState(longint root_id, longint nb_type_id, longint min_state, longint max_state);

  /**
   * Handle signal from FixedPairList that new bond was created.
   */
  void onTupleAdded(longint pid1, longint pid2);
  void onTupleRemoved(longint pid1, longint pid2);

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
  void PrintResTopology();
  void PrintResidues();

  /**
   * Get neighbour list.
   */
  python::list getNeighbourLists();

  static void registerPython();

  /// Defines map that stores res_id -> particle_id
  typedef std::set<longint> PSet;
  typedef std::map<longint, shared_ptr<PSet> > ResParticleIds;

 private:
  typedef std::pair<longint, std::pair<longint, longint> > Triplets;
  typedef std::pair<longint, std::pair<longint, std::pair<longint, longint> > > Quadruplets;
  typedef std::vector<std::pair<longint, longint> > EdgesVector;
  typedef std::vector<std::pair<Particle*, shared_ptr<ParticleProperties> > > NewLocalParticleProperties;
  typedef std::set<std::pair<longint, longint> > SetPairs;
  typedef std::map<longint, longint> MapPairs;
  typedef std::set<longint> SetPids;

  /**
   * Process removing of the bond.
   */
  bool removeBond(longint pid1, longint pid2);

  /**
   * Process new edge in topology.
   */
  void newEdge(longint pid1, longint pid2);

  /**
   * Process removing edge from topology.
   */
  bool deleteEdge(longint pid1, longint pid2);

  /**
   * Update registered FixedTripleList with new entries.
   */
  void defineAngles(std::set<Triplets> &triplets);
  /**
   * Update registered FixedQuadrupleLists with new entries.
   */
  void defineDihedrals(std::set<Quadruplets> &quadruplets);
  void define14tuples(std::set<Quadruplets> &quadruplets);
  /**
   * Update registered FixedTripleList with new entries.
   */
  void undefineAngles(std::set<Triplets> &triplets);
  /**
   * Update registered FixedQuadrupleLists with new entries.
   */
  void undefineDihedrals(std::set<Quadruplets> &quadruplets);
  void undefine14tuples(std::set<Quadruplets> &quadruplets);

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
   * BFS on graph, looking for connected components to update res_id after edge is removed.
   */

  /**
   * Look for nodes at certain distance with BFS algorithm.
   *
   * @param root The root node.
   * @retval The list of node ids.
   */
  std::vector<longint> getNodesAtDistances(longint root);

  void updateParticlePropertiesAtDistance(int id, int distance);

  /**
   * Connecting/Disconnecting to signals.
   */
  void connect();
  void disconnect();

  shared_ptr<System> system_;

  boost::signals2::connection aftIntV2_, aftCalcF_;

  // Mapping for tuples, triplets and quadruplets.
  typedef boost::unordered_map<longint,
      boost::unordered_map<longint, shared_ptr<FixedPairList> > > TupleMap;
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

  bool update_angles_dihedrals;

  std::vector<shared_ptr<FixedPairList> > tuples_;
  std::vector<shared_ptr<FixedPairList> > tuples14_;
  std::vector<shared_ptr<FixedTripleList> > triples_;
  std::vector<shared_ptr<FixedQuadrupleList> > quadruples_;

  TupleMap tupleMap_;
  TupleMap tuple14Map_;
  TripleMap tripleMap_;
  QuadrupleMap quadrupleMap_;

  EdgesVector newEdges_;
  EdgesVector removedEdges_;

  // Residues data.
  std::map<longint, longint> pid_rid;  // particle_id -> res_id;
  GraphMap *residues_;
  void newResEdge(longint first, longint second);

  /** Adjacent list. */
  GraphMap *graph_;
  GraphMap *res_graph_;

  /** Data for DFS */
  longint max_nb_distance_;
  std::set<longint> nb_distances_;
  std::map<longint, std::map<longint, shared_ptr<ParticleProperties> > > distance_type_pp_;
  std::vector<longint> nb_distance_particles_;  //<! Stores the pairs distance; particle_id

  /** Data for bond remove. */
  typedef boost::unordered_map<longint, boost::unordered_set<std::pair<longint, longint> > > DistanceEdges;
  /**
   * Remove edges at distance from the pid (root node).
   *
   * @param pid The root pid.
   * @param edges_to_remove The list of edges to remove.
   */
  void removeNeighbourEdges(size_t pid, SetPairs &edges_to_remove);
  std::set<longint> nb_bond_distances_;
  longint max_bond_nb_distance_;
  std::vector<longint> nb_edges_root_to_remove_;  //<! Stores the pairs: distance; particle_id1, particle_id2
  boost::unordered_map<longint, DistanceEdges> edges_type_distance_pair_types_;

  std::map<longint, shared_ptr<ParticleProperties> > new_type_pp_;
  std::vector<longint> new_local_particle_properties_;
  bool updateParticleProperties(longint pid);

  bool is_dirty_;  ///<! If true then exchangeData will run.

  /** Logger */
  static LOG4ESPP_DECL_LOGGER(theLogger);

  /** Timers. */
  esutil::WallTimer wallTimer;  //!< used for timing

  real timeExchangeData;
  real timeGenerateAnglesDihedrals;
  real timeUpdateNeighbourProperty;
  real timeIsResidueConnected;

  void resetTimers() {
    timeExchangeData = 0.0;
    timeGenerateAnglesDihedrals = 0.0;
    timeUpdateNeighbourProperty = 0.0;
    timeIsResidueConnected = 0.0;
  }

  python::list getTimers();
};

}  // end namespace integrator
}  // end namespace espressoppp
#endif //_INTEGRATOR_TOPOLOGYMANAGER_H