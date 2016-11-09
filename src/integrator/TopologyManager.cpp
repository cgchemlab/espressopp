/*
  Copyright (C) 2015-2016
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

#include "TopologyManager.hpp"

#include <queue>

#include "boost/format.hpp"
#include "storage/Storage.hpp"
#include "iterator/CellListIterator.hpp"
#include "boost/serialization/map.hpp"
#include "boost/serialization/set.hpp"
#include "boost/serialization/shared_ptr.hpp"


namespace espressopp {
namespace integrator {


LOG4ESPP_LOGGER(TopologyManager::theLogger, "TopologyManager");

using namespace espressopp::iterator;

TopologyManager::TopologyManager(shared_ptr<System> system) :
    Extension(system), system_(system) {
  LOG4ESPP_INFO(theLogger, "TopologyManager");
  type = Extension::all;
  graph_ = new GraphMap();
  res_graph_ = new GraphMap();

  residues_ = new GraphMap();

  update_angles_dihedrals = false;

  wallTimer.reset();
}

TopologyManager::~TopologyManager() {
  disconnect();
  // Clean graph data structure
  for (GraphMap::iterator it = graph_->begin(); it != graph_->end(); ++it) {
    if (it->second)
      delete it->second;
  }
  delete graph_;

  for (GraphMap::iterator it = res_graph_->begin(); it != res_graph_->end(); ++it) {
    if (it->second)
      delete it->second;
  }
  delete res_graph_;

  for (GraphMap::iterator it = residues_->begin(); it != residues_->end(); ++it) {
    if (it->second)
      delete it->second;
  }
  delete residues_;
}

void TopologyManager::connect() {
  aftIntV2_ = integrator->aftIntV2.connect(boost::bind(&TopologyManager::exchangeData, this));
}

void TopologyManager::disconnect() {
  aftIntV2_.disconnect();
}

void TopologyManager::observeTuple(shared_ptr<FixedPairList> fpl) {
  fpl->onTupleAdded.connect(
      boost::bind(&TopologyManager::onTupleAdded, this, _1, _2));
  fpl->onTupleRemoved.connect(
      boost::bind(&TopologyManager::onTupleRemoved, this, _1, _2));
  tuples_.push_back(fpl);
}

/** Registers methods, those FixedList are only to updated and no to take data. */
void TopologyManager::registerTuple(
    shared_ptr<FixedPairList> fpl, longint type1, longint type2, longint level) {
  tuples_.push_back(fpl);
  tupleMap_[type1][type2] = fpl;
  tupleMap_[type2][type1] = fpl;
}

void TopologyManager::registerTriple(shared_ptr<FixedTripleList> ftl,
                                     longint type1, longint type2, longint type3) {
  tripleMap_[type1][type2][type3] = ftl;
  tripleMap_[type3][type2][type1] = ftl;
  triples_.push_back(ftl);
  update_angles_dihedrals = true;
}

void TopologyManager::registerQuadruple(shared_ptr<FixedQuadrupleList> fql, longint type1,
                                        longint type2, longint type3, longint type4) {
  quadrupleMap_[type1][type2][type3][type4] = fql;
  quadrupleMap_[type4][type3][type2][type1] = fql;
  quadruples_.push_back(fql);
  update_angles_dihedrals = true;
}

void TopologyManager::InitializeTopology() {
  // Collect locally the list of edges by iterating over registered tuple lists with bonds.
  EdgesVector edges;
  EdgesVector res_ids;
  EdgesVector output;
  std::map<longint, std::set<int> > tmp_resGraph;
  std::vector<std::pair<longint, longint> > local_resid;

  for (std::vector< shared_ptr<FixedPairList> >::iterator it = tuples_.begin(); it != tuples_.end(); ++it) {
    for (FixedPairList::PairList::Iterator pit(**it); pit.isValid(); ++pit) {
      Particle &p1 = *pit->first;
      Particle &p2 = *pit->second;
      edges.push_back(std::make_pair(p1.id(), p2.id()));
      // Manage bond between residues.
      bool foundResBond = false;
      if (tmp_resGraph.count(p1.res_id()) > 0)
        if (tmp_resGraph.at(p1.res_id()).count(p2.res_id()) > 0)
          foundResBond = true;
      else if (tmp_resGraph.count(p2.res_id()) > 0)
          if (tmp_resGraph.at(p2.res_id()).count(p1.res_id()) > 0)
            foundResBond = true;

      if (!foundResBond) {
        res_ids.push_back(std::make_pair(p1.res_id(), p2.res_id()));
        tmp_resGraph[p1.res_id()].insert(p2.res_id());
        tmp_resGraph[p2.res_id()].insert(p1.res_id());
      }
      //local_resid.push_back(std::make_pair(p1.id(), p1.res_id()));
      //local_resid.push_back(std::make_pair(p2.id(), p2.res_id()));
    }
  }
  // Make global map of pid->res_id
  CellList cells = system_->storage->getRealCells();
  for(CellListIterator cit(cells); !cit.isDone(); ++cit) {
    local_resid.push_back(std::make_pair(cit->id(), cit->res_id()));
  }

  output.push_back(std::make_pair(edges.size(), res_ids.size()));
  output.push_back(std::make_pair(local_resid.size(), 0));
  output.insert(output.end(), edges.begin(), edges.end());
  output.insert(output.end(), res_ids.begin(), res_ids.end());
  output.insert(output.end(), local_resid.begin(), local_resid.end());

  LOG4ESPP_DEBUG(theLogger, "Scatter " << output.size());
  // Scatter edges lists all over all nodes. This is costful operation but
  // it is simpler than moving part of graphs all around.
  std::vector<EdgesVector> global_output;
  mpi::all_gather(*(system_->comm), output, global_output);

  // Build a graph. The same on every CPU.
  for (std::vector<EdgesVector>::iterator itv = global_output.begin();
       itv != global_output.end(); ++itv) {
    EdgesVector::iterator it = itv->begin();
    // collects the size of data structures.
    longint edges_size = it->first;
    longint resid_size = it->second;
    it++;
    longint local_resid_size = it->first;
    it++;

    // Create the edge list between atoms.
    for (longint i = 0; i < edges_size; ++it, i++) {
      newEdge(it->first, it->second);
    }

    // Create the edge list between residues.
    for (longint i = 0; i < resid_size; ++it, i++) {
      newResEdge(it->first, it->second);
    }

    // Create the mapping particle_id->residue_id and residue_id->particle_list
    for (longint i = 0; i < local_resid_size; ++it, i++) {
      longint pid = it->first;
      longint rid = it->second;
      pid_rid[pid] = rid;
      //pid_rid.insert(std::make_pair(pid, rid));  // particle_id -> residue_id;
      if (residues_->count(rid) == 0)
        residues_->insert(std::make_pair(rid, new std::set<longint>()));
      residues_->at(rid)->insert(pid);
    }
  }
}

void TopologyManager::PrintTopology() {
  for (GraphMap::iterator it = graph_->begin(); it != graph_->end(); ++it) {
    if (it->second != NULL) {
      std::cout << it->first << ": ";
      for (std::set<int>::iterator itv = it->second->begin(); itv != it->second->end(); ++itv) {
        std::cout << *itv << " ";
      }
      std::cout << std::endl;
    }
  }
}

void TopologyManager::PrintResTopology() {
  for (GraphMap::iterator it = res_graph_->begin(); it != res_graph_->end(); ++it) {
    if (it->second != NULL) {
      std::cout << it->first << ": ";
      for (std::set<int>::iterator itv = it->second->begin(); itv != it->second->end(); ++itv) {
        std::cout << *itv << " ";
      }
      std::cout << std::endl;
    }
  }
}

void TopologyManager::PrintResidues() {
  for (GraphMap::iterator it = residues_->begin(); it != residues_->end(); ++it) {
    if (it->second != NULL) {
      std::cout << it->first << ": ";
      for (std::set<int>::iterator itv = it->second->begin(); itv != it->second->end(); ++itv) {
        std::cout << *itv << " ";
      }
      std::cout << std::endl;
    }
  }

  std::cout << "Map PID->RID" << std::endl;
  for (std::map<longint, longint>::iterator it = pid_rid.begin(); it != pid_rid.end(); ++it) {
    std::cout << it->first << ": " << it->second << std::endl;
  }

}

python::list TopologyManager::getNeighbourLists() {
  python::list nodes;
  for (GraphMap::iterator it = graph_->begin(); it != graph_->end(); ++it) {
    if (it->second != NULL) {
      python::list neighbours;
      for (std::set<int>::iterator itv = it->second->begin(); itv != it->second->end(); ++itv) {
        neighbours.append(*itv);
      }
      nodes.append(python::make_tuple(it->first, neighbours));
    }
  }
  return nodes;
}


void TopologyManager::onTupleAdded(longint pid1, longint pid2) {
  LOG4ESPP_DEBUG(theLogger, "onTupleAdded pid1=" << pid1 << " pid2=" << pid2);
  newEdges_.push_back(std::make_pair(pid1, pid2));
}


void TopologyManager::newEdge(longint pid1, longint pid2) {
  /// Updates graph.
  if (graph_->count(pid1) == 0)
    graph_->insert(std::make_pair(pid1, new std::set<longint>()));
  if (graph_->count(pid2) == 0)
    graph_->insert(std::make_pair(pid2, new std::set<longint>()));
  graph_->at(pid1)->insert(pid2);
  graph_->at(pid2)->insert(pid1);

  // newResidueEdge
  newResEdge(pid_rid[pid1], pid_rid[pid2]);

  // Generate angles, dihedrals, based on updated graph.
  wallTimer.startMeasure();
  if (update_angles_dihedrals) {
    std::set<Quadruplets> *quadruplets = new std::set<Quadruplets>();
    std::set<Triplets> *triplets = new std::set<Triplets>();
    generateAnglesDihedrals(pid1, pid2, *quadruplets, *triplets);

    defineAngles(*triplets);
    defineDihedrals(*quadruplets);
    define14tuples(*quadruplets);
  }
  timeGenerateAnglesDihedrals += wallTimer.stopMeasure();
}

void TopologyManager::newResEdge(longint rpid1, longint rpid2) {
  LOG4ESPP_DEBUG(theLogger, "newResEdge rpid1=" << rpid1 << " rpid2=" << rpid2);
  /// Updates residue graph.
  if (res_graph_->count(rpid1) == 0)
    res_graph_->insert(std::make_pair(rpid1, new std::set<longint>()));
  if (res_graph_->count(rpid2) == 0)
    res_graph_->insert(std::make_pair(rpid2, new std::set<longint>()));
  res_graph_->at(rpid1)->insert(rpid2);
  res_graph_->at(rpid2)->insert(rpid1);
}

void TopologyManager::onTupleRemoved(longint pid1, longint pid2) {

  removedEdges_.push_back(std::make_pair(pid1, pid2));
}

void TopologyManager::deleteEdge(longint pid1, longint pid2) {
  removeBond(pid1, pid2);  // remove bond from fpl

  if (graph_->count(pid1) > 0)
    graph_->at(pid1)->erase(pid2);
  if (graph_->count(pid2) > 0)
    graph_->at(pid2)->erase(pid1);

  // If edge removed, check if there is still edge between residues.
  longint rid1 = pid_rid[pid1];
  longint rid2 = pid_rid[pid2];
  // Get list of particles in given residues.
  std::set<longint> *Pset1;
  if (residues_->find(rid1) != residues_->end())
    Pset1 = residues_->at(rid1);
  else
    throw std::runtime_error((const std::string &) (boost::format("Pset1 residue id %d not found") % rid1));
  std::set<longint> *Pset2;
  if (residues_->find(rid2) != residues_->end())
    Pset2 = residues_->at(rid2);
  else
    throw std::runtime_error((const std::string &) (boost::format("Pset2 residue id %d not found") % rid2));

  // Scan through the bonds that can be created between two residues and check if bond exitst.
  // if not then remove bond between residues.
  bool hasBond = false;
  for (std::set<longint>::iterator it1 = Pset1->begin(); !hasBond && it1 != Pset1->end(); ++it1) {
    for (std::set<longint>::iterator it2 = Pset2->begin(); !hasBond && it2 != Pset2->end(); ++it2) {
      longint ppid1 = *it1;
      longint ppid2 = *it2;
      if (graph_->count(ppid1) == 1 && graph_->at(ppid1)->count(ppid2) == 1)
        hasBond = true;
    }
  }
  if (!hasBond) {
    res_graph_->at(rid1)->erase(rid2);
    res_graph_->at(rid2)->erase(rid1);
  }
}

void TopologyManager::removeBond(longint pid1, longint pid2) {
  LOG4ESPP_DEBUG(theLogger, "Removed bond; removing edge: " << pid1 << "-" << pid2);

  Particle *p1 = system_->storage->lookupLocalParticle(pid1);
  Particle *p2 = system_->storage->lookupLocalParticle(pid2);
  if (p1->ghost() && !p2->ghost()) {
    std::swap(pid1, pid2);
  } else if (p1->ghost() && p2->ghost()) {
    return;
  }

  longint t1 = p1->type();
  longint t2 = p2->type();
  // Update fpl and remove bond.
  shared_ptr<FixedPairList> fpl = tupleMap_[t1][t2];
  if (!fpl)
    fpl = tupleMap_[t2][t1];

  if (fpl) {
    fpl->remove(pid1, pid2);
  } else {
    throw std::runtime_error((const std::string &) (boost::format("Tuple for pair type %d-%d not found") % t1 % t2));
  }

  // Generate list of angles/dihedrals to remove, based on the graph.
  wallTimer.startMeasure();
  if (update_angles_dihedrals) {
    std::set <Quadruplets> *quadruplets = new std::set<Quadruplets>();
    std::set <Triplets> *triplets = new std::set<Triplets>();
    generateAnglesDihedrals(pid1, pid2, *quadruplets, *triplets);

    // Update angles.
    undefineAngles(*triplets);
    undefineDihedrals(*quadruplets);
    undefine14tuples(*quadruplets);
    for (std::vector<shared_ptr<FixedTripleList> >::iterator it = triples_.begin();
        it != triples_.end(); ++it) {
      (*it)->updateParticlesStorage();
    }
    for (std::vector<shared_ptr<FixedQuadrupleList> >::iterator it = quadruples_.begin();
        it != quadruples_.end(); ++it) {
      (*it)->updateParticlesStorage();
    }
  }
  timeGenerateAnglesDihedrals += wallTimer.stopMeasure();
}

/**
 * Invoke by signal from integrator aftIntF (after all other part)
 * Exchange data between nodes and perform operations.
 */
void TopologyManager::exchangeData() {
  LOG4ESPP_DEBUG(theLogger, "entering exchangeData");
  wallTimer.startMeasure();

  // Collect all message from other CPUs. Both for res_id and new graph edges.
  typedef std::vector<std::vector<longint> > GlobalMerge;
  GlobalMerge global_merge_sets;

  // Pack data
  std::vector<longint> output;
  output.push_back(nb_edges_root_to_remove_.size());
  output.push_back(nb_distance_particles_.size() / 2);  // vector of particles to updates.
  output.push_back(newEdges_.size());  // vector of new edges.
  output.push_back(removedEdges_.size());  // vector of edges to remove.

  output.insert(output.end(), nb_edges_root_to_remove_.begin(), nb_edges_root_to_remove_.end());
  output.insert(output.end(), nb_distance_particles_.begin(), nb_distance_particles_.end());

  for (std::vector<std::pair<longint, longint> >::iterator it = newEdges_.begin();
       it != newEdges_.end(); ++it) {
    output.push_back(it->first);
    output.push_back(it->second);
  }

  for (std::vector<std::pair<longint, longint> >::iterator it = removedEdges_.begin();
       it != removedEdges_.end(); ++it) {
    output.push_back(it->first);
    output.push_back(it->second);
  }

  // Send and gather data from all nodes.
  mpi::all_gather(*(system_->comm), output, global_merge_sets);

  // Merge data from other nodes.
  longint f1, f2, nb_edges_root_to_remove_size;
  longint new_edge_size, remove_edge_size, nb_distance_particles_size;
  for (GlobalMerge::iterator gms = global_merge_sets.begin(); gms != global_merge_sets.end(); gms++) {
    for (std::vector<longint>::iterator itm = gms->begin(); itm != gms->end();) {
      nb_edges_root_to_remove_size = *(itm++);
      nb_distance_particles_size = *(itm++);
      new_edge_size = *(itm++);
      remove_edge_size = *(itm++);

      for (int i = 0; i < nb_edges_root_to_remove_size; i++) {
        int particle_id = *(itm++);
        removeNeighbourEdges(particle_id);
      }

      for (int i = 0; i < nb_distance_particles_size; i++) {
        int distance = *(itm++);
        int particle_id = *(itm++);
        updateParticlePropertiesAtDistance(particle_id, distance);
      }
      for (int i = 0; i < new_edge_size; i++) {
        f1 = *(itm++);
        f2 = *(itm++);
        newEdge(f1, f2);
      }
      for (int i = 0; i < remove_edge_size; i++) {
        f1 = *(itm++);
        f2 = *(itm++);
        deleteEdge(f1, f2);
      }
    }
  }
  newEdges_.clear();
  removedEdges_.clear();
  nb_distance_particles_.clear();

  timeExchangeData += wallTimer.stopMeasure();
  LOG4ESPP_DEBUG(theLogger, "leaving exchangeData");
}

void TopologyManager::defineAngles(std::set<Triplets> &triplets) {
  LOG4ESPP_DEBUG(theLogger, "entering update angles");
  longint t1, t2, t3;
  shared_ptr<FixedTripleList> ftl;
  for (std::set<Triplets>::iterator it = triplets.begin(); it != triplets.end(); ++it) {
    Particle *p1 = system_->storage->lookupLocalParticle(it->first);
    Particle *p2 = system_->storage->lookupLocalParticle(it->second.first);
    Particle *p3 = system_->storage->lookupLocalParticle(it->second.second);
    if (p1 && p2 && p3) {
      t1 = p1->type();
      t2 = p2->type();
      t3 = p3->type();
      // Look for fixed triple list which should be updated.
      ftl = tripleMap_[t1][t2][t3];
      if (!ftl)
        ftl = tripleMap_[t3][t2][t1];
      if (ftl) {
        LOG4ESPP_DEBUG(theLogger, "Found tuple for: " << t1 << "-" << t2 << "-" << t3);
        bool ret = ftl->iadd(p1->id(), p2->id(), p3->id());
        if (!ret)
          ret = ftl->iadd(p3->id(), p2->id(), p1->id());
        if (ret) LOG4ESPP_DEBUG(theLogger,
                                "Defined new angle: " << it->first << "-" << it->second.first << "-"
                                    << it->second.second);
      }
    }
  }
  LOG4ESPP_DEBUG(theLogger, "leaving update angles");
}

void TopologyManager::defineDihedrals(std::set<Quadruplets> &quadruplets) {
  LOG4ESPP_DEBUG(theLogger, "entering update dihedrals");
  longint t1, t2, t3, t4;
  shared_ptr<FixedQuadrupleList> fql;
  for (std::set<Quadruplets>::iterator it = quadruplets.begin(); it != quadruplets.end(); ++it) {
    Particle *p1 = system_->storage->lookupLocalParticle(it->first);
    Particle *p2 = system_->storage->lookupLocalParticle(it->second.first);
    Particle *p3 = system_->storage->lookupLocalParticle(it->second.second.first);
    Particle *p4 = system_->storage->lookupLocalParticle(it->second.second.second);
    if (p1 && p2 && p3 && p4) {
      t1 = p1->type();
      t2 = p2->type();
      t3 = p3->type();
      t4 = p4->type();
      // Look for fixed triple list which should be updated.
      fql = quadrupleMap_[t1][t2][t3][t4];
      if (!fql)
        fql = quadrupleMap_[t4][t3][t2][t1];
      if (fql) {
        LOG4ESPP_DEBUG(theLogger, "Found tuple for: " << t1 << "-" << t2 << "-" << t3 << "-" << t4);
        bool ret = fql->iadd(p1->id(), p2->id(), p3->id(), p4->id());
        if (!ret)
          ret = fql->iadd(p4->id(), p3->id(), p2->id(), p1->id());
        if (ret) LOG4ESPP_DEBUG(theLogger,
                                "Defined new dihedral: " << it->first << "-" << it->second.first
                                    << "-" << it->second.second.first << "-"
                                    << it->second.second.second);
      }
    }
  }
  LOG4ESPP_DEBUG(theLogger, "leaving update dihedrals");
}

void TopologyManager::define14tuples(std::set<Quadruplets> &quadruplets) {
  LOG4ESPP_DEBUG(theLogger, "entering define 1-4 tuples");
  longint t1, t4;
  shared_ptr<FixedPairList> fpl;
  for (std::set<Quadruplets>::iterator it = quadruplets.begin(); it != quadruplets.end(); ++it) {
    Particle *p1 = system_->storage->lookupLocalParticle(it->first);
    Particle *p4 = system_->storage->lookupLocalParticle(it->second.second.second);
    if (p1 && p4) {
      t1 = p1->type();
      t4 = p4->type();
      // Look for fixed triple list which should be updated.
      fpl = tupleMap_[t1][t4];
      if (!fpl)
        fpl = tupleMap_[t4][t1];
      if (fpl) {
        LOG4ESPP_DEBUG(theLogger, "Found tuple for: " << t1 << "-" << t4);
        bool ret = fpl->iadd(p1->id(), p4->id());
        if (!ret)
          ret = fpl->iadd(p4->id(), p1->id());
        if (ret) LOG4ESPP_DEBUG(theLogger,
                                "Defined new 1-4 pair: " << it->first << "-"
                                    << it->second.second.second);
      }
    }
  }
}

void TopologyManager::undefineAngles(std::set<Triplets> &triplets) {
  LOG4ESPP_DEBUG(theLogger, "entering undefine angles");
  longint t1, t2, t3;
  shared_ptr<FixedTripleList> ftl;
  for (std::set<Triplets>::iterator it = triplets.begin(); it != triplets.end(); ++it) {
    Particle *p1 = system_->storage->lookupLocalParticle(it->first);
    Particle *p2 = system_->storage->lookupLocalParticle(it->second.first);
    Particle *p3 = system_->storage->lookupLocalParticle(it->second.second);
    if (p1 && p2 && p3) {
      t1 = p1->type();
      t2 = p2->type();
      t3 = p3->type();
      // Look for fixed triple list which should be updated.
      ftl = tripleMap_[t1][t2][t3];
      if (!ftl)
        ftl = tripleMap_[t3][t2][t1];
      if (ftl) {
        LOG4ESPP_DEBUG(theLogger, "Found tuple for: " << t1 << "-" << t2 << "-" << t3);
        bool ret = ftl->remove(p1->id(), p2->id(), p3->id());
        if (!ret)
          ret = ftl->remove(p3->id(), p2->id(), p1->id());
        if (ret) {
          LOG4ESPP_DEBUG(theLogger,
                         "Remove angle: " << it->first << "-" << it->second.first << "-"
                             << it->second.second);
        }
      }
    }
  }
  LOG4ESPP_DEBUG(theLogger, "leaving undefine angles");
}

void TopologyManager::undefineDihedrals(std::set<Quadruplets> &quadruplets) {
  longint t1, t2, t3, t4;
  shared_ptr<FixedQuadrupleList> fql;
  for (std::set<Quadruplets>::iterator it = quadruplets.begin(); it != quadruplets.end(); ++it) {
    Particle *p1 = system_->storage->lookupLocalParticle(it->first);
    Particle *p2 = system_->storage->lookupLocalParticle(it->second.first);
    Particle *p3 = system_->storage->lookupLocalParticle(it->second.second.first);
    Particle *p4 = system_->storage->lookupLocalParticle(it->second.second.second);
    if (p1 && p2 && p3 && p4) {
      t1 = p1->type();
      t2 = p2->type();
      t3 = p3->type();
      t4 = p4->type();
      // Look for fixed triple list which should be updated.
      fql = quadrupleMap_[t1][t2][t3][t4];
      if (!fql)
        fql = quadrupleMap_[t4][t3][t2][t1];
      if (fql) {
        LOG4ESPP_DEBUG(theLogger, "Found tuple for: " << t1 << "-" << t2 << "-" << t3 << "-" << t4);
        bool ret = fql->remove(p1->id(), p2->id(), p3->id(), p4->id());

        if (!ret)
          ret = fql->remove(p4->id(), p3->id(), p2->id(), p1->id());
        if (ret) {
          LOG4ESPP_DEBUG(theLogger,
                         "Remove dihedral: " << it->first << "-" << it->second.first
                             << "-" << it->second.second.first << "-"
                             << it->second.second.second);
        }
      }
    }
  }
}

void TopologyManager::undefine14tuples(std::set<Quadruplets> &quadruplets) {
  longint t1, t4;
  shared_ptr<FixedPairList> fpl;
  for (std::set<Quadruplets>::iterator it = quadruplets.begin(); it != quadruplets.end(); ++it) {
    Particle *p1 = system_->storage->lookupLocalParticle(it->first);
    Particle *p4 = system_->storage->lookupLocalParticle(it->second.second.second);
    if (p1 && p4) {
      t1 = p1->type();
      t4 = p4->type();
      // Look for fixed triple list which should be updated.
      fpl = tupleMap_[t1][t4];
      if (!fpl)
        fpl = tupleMap_[t4][t1];
      if (fpl) {
        LOG4ESPP_DEBUG(theLogger, "Found tuple for: " << t1 << "-" << t4);
        bool ret = fpl->remove(p1->id(), p4->id());
        if (!ret)
          ret = fpl->remove(p4->id(), p1->id());
        if (ret) LOG4ESPP_DEBUG(theLogger,
                                "Remove dihedral: " << it->first << "-"
                                    << it->second.second.second);
      }
    }
  }
}

void TopologyManager::generateAnglesDihedrals(longint pid1,
                                              longint pid2,
                                              std::set<Quadruplets> &quadruplets,
                                              std::set<Triplets> &triplets) {
  std::set<longint> *nb1, *nb2;
  if (graph_->find(pid1) != graph_->end())
    nb1 = graph_->at(pid1);
  else
    throw std::runtime_error((const std::string &) (boost::format("Node pid1 %d not found") % pid1));
  if (graph_->find(pid2) != graph_->end())
    nb2 = graph_->at(pid2);
  else
    throw std::runtime_error((const std::string &) (boost::format("Node pid2 %d not found") % pid1));
  // Case pid2 pid1 <> <>
  if (nb1) {
    //Iterates over p1 neighbours
    for (std::set<longint>::iterator it = nb1->begin(); it != nb1->end(); ++it) {
      if (*it == pid1 || *it == pid2)
        continue;
      if (triplets.count(std::make_pair(*it, std::make_pair(pid2, pid1))) == 0)
        triplets.insert(std::make_pair(pid2, std::make_pair(pid1, *it)));

      std::set<longint> *nbb1 = graph_->at(*it);
      if (nbb1) {
        for (std::set<longint>::iterator itt = nbb1->begin(); itt != nbb1->end(); ++itt) {
          if (*itt == *it || *itt == pid1 || *itt == pid2)
            continue;
          if (quadruplets.count(std::make_pair(*itt,
                                               std::make_pair(*it, std::make_pair(pid1, pid2))))
              == 0)
            quadruplets.insert(
                std::make_pair(pid2, std::make_pair(pid1, std::make_pair(*it, *itt))));
        }
      }
    }
  }
  // Case pid1 pid2 <> <>
  if (nb2) {
    for (std::set<longint>::iterator it = nb2->begin(); it != nb2->end(); ++it) {
      if (*it == pid1 || *it == pid2)
        continue;
      if (triplets.count(std::make_pair(*it, std::make_pair(pid2, pid1))) == 0)
        triplets.insert(std::make_pair(pid1, std::make_pair(pid2, *it)));

      std::set<longint> *nbb2 = graph_->at(*it);
      if (nbb2) {
        for (std::set<longint>::iterator itt = nbb2->begin(); itt != nbb2->end(); ++itt) {
          if (*itt == *it || *itt == pid1 || *itt == pid2)
            continue;
          if (quadruplets.count(std::make_pair(*itt,
                                               std::make_pair(*it, std::make_pair(pid2, pid1))))
              == 0)
            quadruplets.insert(
                std::make_pair(pid1, std::make_pair(pid2, std::make_pair(*it, *itt))));
        }
      }
    }
  }
  // Case <> pid1 pid2 <>
  if (nb1 && nb2) {
    for (std::set<longint>::iterator it1 = nb1->begin(); it1 != nb1->end(); ++it1) {
      if (*it1 == pid1 || *it1 == pid2)
        continue;
      for (std::set<longint>::iterator it2 = nb2->begin(); it2 != nb2->end(); ++it2) {
        if (*it2 == pid1 || *it2 == pid2 || *it1 == *it2)
          continue;
        if (quadruplets.count(std::make_pair(*it2,
                                             std::make_pair(pid2, std::make_pair(pid1, *it1))))
            == 0)
          quadruplets.insert(
              std::make_pair(*it1, std::make_pair(pid1, std::make_pair(pid2, *it2))));
      }
    }
  }
}

std::vector<longint> TopologyManager::getNodesAtDistances(longint root) {
  std::map<longint, longint> visitedDistance;
  std::queue<longint> Q;
  Q.push(root);
  visitedDistance.insert(std::make_pair(root, 0));

  std::vector<longint> nb_at_distance;

  longint current, node, new_distance;
  while (!Q.empty()) {
    current = Q.front();
    new_distance = visitedDistance[current] + 1;
    try {
      std::set<longint> *adj = graph_->at(current);
      for (std::set<longint>::iterator ia = adj->begin(); ia != adj->end(); ++ia) {
        node = *ia;
        if (visitedDistance.count(node) == 0) {
          if (nb_distances_.count(new_distance) == 1) {
            nb_at_distance.push_back(new_distance);
            nb_at_distance.push_back(node);
          }
          if (new_distance < max_nb_distance_) {
            Q.push(node);
          }
          visitedDistance.insert(std::make_pair(node, new_distance));
        }
      }
    } catch (...) {
      std::cout << "Exception graph_->at(current) = " << current << std::endl;
      throw new std::runtime_error("Exception graph_->at");
    }
    Q.pop();
  }
  return nb_at_distance;
}

void TopologyManager::removeNeighbourEdges(size_t pid) {

  std::map<longint, longint> visitedDistance;
  std::queue<longint> Q;

  Particle *root = system_->storage->lookupLocalParticle(pid);
  if (!root)
    return;
  Q.push(root->id());
  visitedDistance.insert(std::make_pair(root->id(), 0));

  boost::unordered_map<longint, DistanceEdges>::iterator distance_edges = edges_type_distance_pair_types_.find(
      root->type());

  if (distance_edges == edges_type_distance_pair_types_.end())
    return;

  boost::unordered_set<std::pair<longint, longint> > pair_types_at_distance;
  DistanceEdges::iterator pair_types_at_distance_iter_;
  boost::unordered_set<std::pair<longint, longint> > edges_to_remove;

  longint current_node, node, new_distance;
  longint type_p1 = -1;
  longint type_p2 = -1;
  while (!Q.empty()) {
    current_node = Q.front();
    Particle *p1_current = system_->storage->lookupLocalParticle(current_node);
    if (p1_current)
      type_p1 = p1_current->type();
    new_distance = visitedDistance[current_node] + 1;
    pair_types_at_distance_iter_ = distance_edges->second.find(new_distance);
    try {
      std::set<longint> *adj = graph_->at(current_node);
      bool has_pairs_at_distance = false;
      if (pair_types_at_distance_iter_ != distance_edges->second.end()) {
        pair_types_at_distance = pair_types_at_distance_iter_->second;
        has_pairs_at_distance = true;
      }
      for (std::set<longint>::iterator ia = adj->begin(); ia != adj->end(); ++ia) {
        node = *ia;
        if (visitedDistance.count(node) == 0) {
          if (has_pairs_at_distance) {
            Particle *p1_node = system_->storage->lookupLocalParticle(node);
            if (p1_node && p1_current) {
              type_p2 = p1_node->type();
              if (pair_types_at_distance.count(std::make_pair(type_p1, type_p2)) != 0) {
                if (edges_to_remove.count(std::make_pair(node, current_node)) == 0)
                  edges_to_remove.insert(std::make_pair(current_node, node));
              }
            }
          }
          if (new_distance < max_bond_nb_distance_) {
            Q.push(node);
          }
          visitedDistance.insert(std::make_pair(node, new_distance));
        }
      }
    } catch (...) {
      std::cout << "Exception graph_->at(current) = " << current_node << std::endl;
      throw new std::runtime_error("Exception graph_->at");
    }
    Q.pop();
  }

  if (edges_to_remove.size() > 0) {
    LOG4ESPP_DEBUG(theLogger, "edges to remove: " << edges_to_remove.size());
    for (boost::unordered_set<std::pair<longint, longint> >::iterator it = edges_to_remove.begin();
         it != edges_to_remove.end(); ++it) {
      deleteEdge(it->first, it->second);
    }
  }
}

void TopologyManager::registerNeighbourPropertyChange(
      longint type_id, shared_ptr<ParticleProperties> pp, longint nb_level) {
  LOG4ESPP_DEBUG(theLogger, "register property change for type_id=" << type_id
      << " at level=" << nb_level);
  max_nb_distance_ = std::max(max_nb_distance_, nb_level);
  nb_distances_.insert(nb_level);
  distance_type_pp_[nb_level][type_id] = pp;
}

void TopologyManager::registerNeighbourBondToRemove(longint type_id,
                                                    longint nb_level,
                                                    longint type_pid1,
                                                    longint type_pid2) {
  max_bond_nb_distance_ = std::max(max_bond_nb_distance_, nb_level);
  edges_type_distance_pair_types_[type_id][nb_level].insert(std::make_pair(type_pid1, type_pid2));
  edges_type_distance_pair_types_[type_id][nb_level].insert(std::make_pair(type_pid2, type_pid1));
}


void TopologyManager::invokeNeighbourPropertyChange(Particle &root) {
  wallTimer.startMeasure();
  std::vector<longint> nb = getNodesAtDistances(root.id());
  LOG4ESPP_DEBUG(theLogger, "inovokeNeighbourPropertyChange from root=" << root.id()
      << " generates=" << nb.size() << " of neighbour particles");
  nb_distance_particles_.insert(nb_distance_particles_.end(), nb.begin(), nb.end());
  timeUpdateNeighbourProperty += wallTimer.stopMeasure();
}

void TopologyManager::invokeNeighbourBondRemove(Particle &root) {
  // Check if this root.type is on the list of possible edges to remove.
  if (edges_type_distance_pair_types_.count(root.type()) == 1) {
    nb_edges_root_to_remove_.insert(nb_edges_root_to_remove_.end(), root.id());
  }
}

void TopologyManager::updateParticlePropertiesAtDistance(int pid, int distance) {
  LOG4ESPP_DEBUG(theLogger, "update particle properties id=" << pid << " at distance=" << distance);
  // We will update both ghost and normal particles as ghost can also take part in reactions.
  Particle *p = system_->storage->lookupLocalParticle(pid);

  if (p) {  // particle exists here.
    longint p_type = p->type();
    if (distance_type_pp_.count(distance) > 0) {
      if (distance_type_pp_[distance].count(p_type) > 0) {
        shared_ptr<ParticleProperties> pp = distance_type_pp_[distance][p_type];
        pp->updateParticleProperties(p);
      }
    }
  }
}

bool TopologyManager::isResiduesConnected(longint rid1, longint rid2) {
  wallTimer.startMeasure();
  bool ret = (res_graph_->count(rid1) == 1 && res_graph_->at(rid1)->count(rid2) == 1);
  timeIsResidueConnected += wallTimer.stopMeasure();
  return ret;
}

bool TopologyManager::isParticleConnected(longint pid1, longint pid2) {
  bool ret = (graph_->count(pid1) == 1 && graph_->at(pid1)->count(pid2) == 1);
  return ret;
}


void TopologyManager::registerPython() {
  using namespace espressopp::python;

  boost::python::implicitly_convertible<shared_ptr<FixedPairListLambda>, shared_ptr<FixedPairList> >();
  boost::python::implicitly_convertible<shared_ptr<FixedTripleListLambda>, shared_ptr<FixedTripleList> >();
  boost::python::implicitly_convertible<shared_ptr<FixedQuadrupleListLambda>, shared_ptr<FixedQuadrupleList> >();

  class_<TopologyManager, shared_ptr<TopologyManager>, bases<Extension> >
      ("integrator_TopologyManager", init<shared_ptr<System> >())
      .def("connect", &TopologyManager::connect)
      .def("disconnect", &TopologyManager::disconnect)
      .def("observe_tuple", &TopologyManager::observeTuple)
      .def("register_tuple", &TopologyManager::registerTuple)
      .def("register_triple", &TopologyManager::registerTriple)
      .def("register_quadruple", &TopologyManager::registerQuadruple)
      .def("initialize", &TopologyManager::InitializeTopology)
      .def("exchange_data", &TopologyManager::exchangeData)
      .def("print_topology", &TopologyManager::PrintTopology)
      .def("print_res_topology", &TopologyManager::PrintResTopology)
      .def("print_residues", &TopologyManager::PrintResidues)
      .def("get_neighbour_lists", &TopologyManager::getNeighbourLists)
      .def("get_timers", &TopologyManager::getTimers)
      .def("is_residue_connected", &TopologyManager::isResiduesConnected)
      .def("is_particle_connected", &TopologyManager::isParticleConnected)
      ;
}
python::list TopologyManager::getTimers() {
  python::list ret;
  ret.append(python::make_tuple("timeExchangeData", timeExchangeData));
  ret.append(python::make_tuple("timeGenerateAnglesDihedrals", timeGenerateAnglesDihedrals));
  ret.append(python::make_tuple("timeUpdateNeighbourProperty", timeUpdateNeighbourProperty));
  ret.append(python::make_tuple("timeIsResidueConnected", timeIsResidueConnected));

  return ret;
}


}  // end namespace integrator
}  // end namespace espressoppp
