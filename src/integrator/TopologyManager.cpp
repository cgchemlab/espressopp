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

#include <algorithm>
#include <queue>
#include <utility>

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
  update_angles_dihedrals = false;
}

TopologyManager::~TopologyManager() {
  disconnect();
  // Clean graph data structure
  for (GraphMap::iterator it = graph_->begin(); it != graph_->end(); ++it) {
    if (it->second)
      delete it->second;
  }
  delete graph_;
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

void TopologyManager::registerTuple(
    shared_ptr<FixedPairList> fpl, longint type1, longint type2, longint level) {
  tuples_.push_back(fpl);
  LOG4ESPP_ERROR(theLogger, "Currently, only calculated 1-4 interaction.");
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
  for (std::vector< shared_ptr<FixedPairList> >::iterator it = tuples_.begin();
       it != tuples_.end(); ++it) {
    for (FixedPairList::PairList::Iterator pit(**it); pit.isValid(); ++pit) {
      Particle &p1 = *pit->first;
      Particle &p2 = *pit->second;
      edges.push_back(std::make_pair(p1.id(), p2.id()));
    }
  }
  LOG4ESPP_DEBUG(theLogger, "Scatter " << edges.size());
  // Scatter edges lists all over all nodes. This is costful operation but
  // it is simpler than moving part of graphs all around.
  std::vector<EdgesVector> global_edges;
  mpi::all_gather(*(system_->comm), edges, global_edges);

  // Build a graph. The same on every CPU.
  for (std::vector<EdgesVector>::iterator itv = global_edges.begin();
       itv != global_edges.end(); ++itv) {
    for (EdgesVector::iterator it = itv->begin(); it != itv->end(); ++it) {
      newEdge(it->first, it->second);
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
  Particle *p1 = system_->storage->lookupRealParticle(pid1);
  Particle *p2 = system_->storage->lookupLocalParticle(pid2);
  if (!p1 || !p2)
    return;
  // Update res_id, if bond joins two molecules with different res_id then merge the sets.
  // p2.res_id() <- p1.res_id()
  if (p1->res_id() != p2->res_id()) {
    merge_sets_.push_back(std::make_pair(p1->res_id(), p2->res_id()));
  }
  newEdges_.push_back(std::make_pair(pid1, pid2));
}


void TopologyManager::newEdge(longint pid1, longint pid2) {
  if (graph_->count(pid1) == 0)
    graph_->insert(std::make_pair(pid1, new std::set<longint>()));
  if (graph_->count(pid2) == 0)
    graph_->insert(std::make_pair(pid2, new std::set<longint>()));
  graph_->at(pid1)->insert(pid2);
  graph_->at(pid2)->insert(pid1);

  // Generate angles, dihedrals, based on updated graph.
  if (update_angles_dihedrals) {
    std::set<Quadruplets> *quadruplets = new std::set<Quadruplets>();
    std::set<Triplets> *triplets = new std::set<Triplets>();
    generateAnglesDihedrals(pid1, pid2, *quadruplets, *triplets);

    defineAngles(*triplets);
    defineDihedrals(*quadruplets);
    define14tuples(*quadruplets);
  }
}

void TopologyManager::onTupleRemoved(longint pid1, longint pid2) {
  Particle *p1 = system_->storage->lookupRealParticle(pid1);
  Particle *p2 = system_->storage->lookupLocalParticle(pid2);
  if (!p1 || !p2)
    return;
  if (p1->res_id() == p2->res_id()) {
    split_sets_.push_back(std::make_pair(p1->res_id(), std::make_pair(pid1, pid2)));
  }
  removedEdges_.push_back(std::make_pair(pid1, pid2));
}

void TopologyManager::deleteEdge(longint pid1, longint pid2) {
  removeBond(pid1, pid2);

  if (graph_->count(pid1) > 0)
    graph_->at(pid1)->erase(pid2);
  if (graph_->count(pid2) > 0)
    graph_->at(pid2)->erase(pid1);
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

  // Generate list of angles/dihedrals to remove, based on the graph.
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
        if (ret) LOG4ESPP_DEBUG(theLogger,
                                "Remove angle: " << it->first << "-" << it->second.first << "-"
                                    << it->second.second);
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
        if (ret) LOG4ESPP_DEBUG(theLogger,
                                "Remove dihedral: " << it->first << "-" << it->second.first
                                    << "-" << it->second.second.first << "-"
                                    << it->second.second.second);
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
  std::set<longint> *nb1 = graph_->at(pid1);
  std::set<longint> *nb2 = graph_->at(pid2);
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

void TopologyManager::exchangeData() {
  LOG4ESPP_DEBUG(theLogger, "entering exchangeData");
  // Collect all message from other CPUs. Both for res_id and new graph edges.
  typedef std::vector<std::vector<longint> > GlobalMerge;
  GlobalMerge global_merge_sets;

  // Pack data
  // Format: vector of pairs.
  // 0: size of merge_sets, size of new_edges;
  // 1: size of split_sets, size of removed_edges;
  // 2..size_of_merge_sets
  // size_of_merge_sets+2...size_of_merge_sets+size_of_new_edges
  //
  std::vector<longint> output;
  output.push_back(merge_sets_.size());  // vector of sets of particles to merge.
  output.push_back(nb_distance_particles_.size() / 2);  // vector of particles to updates.
  output.push_back(newEdges_.size());  // vector of new edges.
  output.push_back(removedEdges_.size());  // vector of edges to remove.
  output.push_back(split_sets_.size());  // vector of sets to split.

  for (std::vector<std::pair<longint, longint> >::iterator it = merge_sets_.begin();
      it != merge_sets_.end(); ++it) {
    output.push_back(it->first);
    output.push_back(it->second);
  }

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

  for (std::vector<std::pair<longint, std::pair<longint, longint> > >::iterator it = split_sets_.begin();
       it != split_sets_.end(); ++it) {
    output.push_back(it->first);
    output.push_back(it->second.first);
    output.push_back(it->second.second);
  }

  // Send and gather data from all nodes.
  mpi::all_gather(*(system_->comm), output, global_merge_sets);

  // Merge data from other nodes.
  longint f1, f2, f3, merge_set_size, new_edge_size, split_set_size, remove_edge_size, nb_distance_particles_size;
  for (GlobalMerge::iterator gms = global_merge_sets.begin(); gms != global_merge_sets.end(); gms++) {
    for (std::vector<longint>::iterator itm = gms->begin(); itm != gms->end();) {
      merge_set_size = *(itm++);
      nb_distance_particles_size = *(itm++);
      new_edge_size = *(itm++);
      split_set_size = *(itm++);
      remove_edge_size = *(itm++);
      for (int i = 0; i < merge_set_size; i++) {
        f1 = *(itm++);
        f2 = *(itm++);
        mergeResIdSets(f1, f2);
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
      for (int i = 0; i < split_set_size; i++) {
        f1 = *(itm++);
        f2 = *(itm++);
        f3 = *(itm++);
        splitResIdSets(f1, f2, f3);
      }
    }
  }
  merge_sets_.clear();
  newEdges_.clear();
  removedEdges_.clear();
  split_sets_.clear();
  nb_distance_particles_.clear();
  LOG4ESPP_DEBUG(theLogger, "leaving exchangeData");
}

void TopologyManager::mergeResIdSets(longint res_id_a, longint res_id_b) {
  LOG4ESPP_DEBUG(theLogger, "merge sets A:" << res_id_a << " with B:" << res_id_b);
  if (res_id_a == res_id_b)
    return;

  if (res_particle_ids_[res_id_b]->begin() == res_particle_ids_[res_id_b]->end())
    return;

  shared_ptr<PSet> setA = res_particle_ids_[res_id_a];

  shared_ptr<PSet> setB = res_particle_ids_[res_id_b];
  assert (setB != NULL);
  // Merge two sets.
  for (PSet::iterator it = setB->begin(); it != setB->end(); ++it) {
    setA->insert(*it);
    Particle *p = system_->storage->lookupLocalParticle(*it);
    if (p) {
      p->setResId(res_id_a);
    }
  }
  res_particle_ids_[res_id_b] = res_particle_ids_[res_id_a];
  LOG4ESPP_DEBUG(theLogger, "leaving merge sets A:" << res_id_a << " with B:" << res_id_b);
  //delete setB;  // free memory as this set was copied to res_id_a.
}

void TopologyManager::splitResIdSets(longint res_id, longint pid1, longint pid2) {
  LOG4ESPP_DEBUG(theLogger, "spliting set " << res_id << " pid1=" << pid1 << " pid2=" << pid2);

  shared_ptr<PSet> setA = res_particle_ids_[res_id];
  longint max_res_id = (--res_particle_ids_.end())->first;

  PSet set_1;
  PSet set_2;

  std::cout << "splitResIdSets not implemented!!!" << std::endl;

}

void TopologyManager::Rebuild() {
  LOG4ESPP_DEBUG(theLogger, "entering Rebuild");
  for (ResParticleIds::iterator it = res_particle_ids_.begin();
       it != res_particle_ids_.end(); ++it) {
    delete &(*it->second);
  }
  res_particle_ids_.clear();
  // Collect local particle ids.
  CellList cl = system_->storage->getRealCells();
  for (CellListIterator it(cl); it.isValid(); ++it) {
    Particle &p = *it;
    if (!res_particle_ids_[p.res_id()])
      res_particle_ids_[p.res_id()] = make_shared<PSet>();
    res_particle_ids_[p.res_id()]->insert(p.id());
  }
  // Sync among CPUs.
  std::vector<ResParticleIds> global_res_particle_ids;
  mpi::all_gather(*(system_->comm), res_particle_ids_, global_res_particle_ids);
  // Update local storage.
  for (std::vector<ResParticleIds>::iterator it = global_res_particle_ids.begin();
       it != global_res_particle_ids.end(); ++it) {
    for (ResParticleIds::iterator itp = it->begin(); itp != it->end(); ++itp) {
      for (PSet::iterator itps = itp->second->begin(); itps != itp->second->end(); ++itps) {
        if (!res_particle_ids_[itp->first])
          res_particle_ids_[itp->first] = make_shared<PSet>();
        res_particle_ids_[itp->first]->insert(*itps);
      }
    }
  }
  LOG4ESPP_DEBUG(theLogger, "leaving Rebuild");
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
    Q.pop();
  }
  return nb_at_distance;
}


void TopologyManager::registerPython() {
  using namespace espressopp::python;

  class_<TopologyManager, shared_ptr<TopologyManager>, bases<Extension> >
      ("integrator_TopologyManager", init<shared_ptr<System> >())
      .def("connect", &TopologyManager::connect)
      .def("disconnect", &TopologyManager::disconnect)
      .def("rebuild", &TopologyManager::Rebuild)
      .def("observe_tuple", &TopologyManager::observeTuple)
      .def("register_tuple", &TopologyManager::registerTuple)
      .def("register_triple", &TopologyManager::registerTriple)
      .def("register_quadruple", &TopologyManager::registerQuadruple)
      .def("initialize", &TopologyManager::InitializeTopology)
      .def("exchange_data", &TopologyManager::exchangeData)
      .def("print_topology", &TopologyManager::PrintTopology)
      .def("get_neighbour_lists", &TopologyManager::getNeighbourLists);
}


void TopologyManager::registerNeighbourPropertyChange(
      longint type_id, shared_ptr<ParticleProperties> pp, longint nb_level) {
  LOG4ESPP_DEBUG(theLogger, "register property change for type_id=" << type_id
      << " at level=" << nb_level);
  max_nb_distance_ = std::max(max_nb_distance_, nb_level);
  nb_distances_.insert(nb_level);
  distance_type_pp_[nb_level][type_id] = pp;
}


void TopologyManager::invokeNeighbourPropertyChange(Particle &root) {
  std::vector<longint> nb = getNodesAtDistances(root.id());
  LOG4ESPP_DEBUG(theLogger, "inokgeNeighbourPropertyChange from root=" << root.id()
      << " generates=" << nb.size() << " of neighbour particles");
  nb_distance_particles_.insert(nb_distance_particles_.end(), nb.begin(), nb.end());
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

        // Change particle type.
        if (pp->type != NULL) {
          p->setType(pp->type);
        }

        // Change mass.
        if (pp->mass != NULL) {
          p->setMass(pp->mass);
        }

        // Change partial charge.
        if (pp->q != NULL) {
          p->setQ(pp->q);
        }
      }
    }
  }
}

}  // end namespace integrator
}  // end namespace espressoppp
