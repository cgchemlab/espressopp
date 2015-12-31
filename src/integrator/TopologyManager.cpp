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

#include "TopologyManager.hpp"
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

void TopologyManager::observeTuple(shared_ptr<FixedPairList> fpl, longint type1, longint type2) {
  fpl->onTupleAdded.connect(
      boost::bind(&TopologyManager::onTupleAdded, this, _1, _2));
  tupleMap_.push_back(fpl);
}

void TopologyManager::observeTriple(shared_ptr<FixedTripleList> ftl,
                                    longint type1, longint type2, longint type3) {
  tripleMap_[type1][type2][type3] = ftl;
  tripleMap_[type3][type2][type1] = ftl;
  update_angles_dihedrals = true;
}

void TopologyManager::observeQuadruple(shared_ptr<FixedQuadrupleList> fql, longint type1,
                                       longint type2, longint type3, longint type4) {
  quadrupleMap_[type1][type2][type3][type4] = fql;
  quadrupleMap_[type4][type3][type2][type1] = fql;
  update_angles_dihedrals = true;
}

void TopologyManager::InitializeTopology() {
  // Collect locally the list of edges by iterating over registered tuple lists with bonds.
  EdgesVector edges;
  for (std::vector<shared_ptr<FixedPairList> >::iterator it = tupleMap_.begin();
       it != tupleMap_.end(); ++it) {
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
  newBond(pid1, pid2);
}

void TopologyManager::newEdge(longint pid1, longint pid2) {
  if (graph_->count(pid1) == 0)
    graph_->insert(std::make_pair(pid1, new std::set<longint>()));
  if (graph_->count(pid2) == 0)
    graph_->insert(std::make_pair(pid2, new std::set<longint>()));
  graph_->at(pid1)->insert(pid2);
  graph_->at(pid2)->insert(pid1);
}

void TopologyManager::newBond(longint pid1, longint pid2) {
  LOG4ESPP_DEBUG(theLogger, "New bond; Adding edge: " << pid1 << "-" << pid2);
  newEdge(pid1, pid2);

  // Because graph is updated.
  newEdges_.push_back(std::make_pair(pid1, pid2));

  // Generate angles and dihedrals.
  if (update_angles_dihedrals) {
    std::set<Quadruplets> *quadruplets = new std::set<Quadruplets>();
    std::set<Triplets> *triplets = new std::set<Triplets>();
    generateAnglesDihedrals(pid1, pid2, *quadruplets, *triplets);

    updateAngles(*triplets);
    updateDihedrals(*quadruplets);
  }
}

void TopologyManager::updateAngles(std::set<Triplets> &triplets) {
  LOG4ESPP_DEBUG(theLogger, "entering update angles");
  longint t1, t2, t3;
  shared_ptr<FixedTripleList> ftl;
  for (std::set<Triplets>::iterator it = triplets.begin(); it != triplets.end(); ++it) {
    Particle *p1 = system_->storage->lookupLocalParticle(it->first);
    Particle *p2 = system_->storage->lookupLocalParticle(it->second.first);
    Particle *p3 = system_->storage->lookupLocalParticle(it->second.second);
    if ((p1 && p2) && (p1 && p3)) {
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

void TopologyManager::updateDihedrals(std::set<Quadruplets> &quadruplets) {
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
  typedef std::vector<std::vector<std::pair<longint, longint> > > GlobaleMergeSets;
  GlobaleMergeSets global_merge_sets;

  // Format: vector of pairs.
  // 0: size of merge_sets, size of new_edges
  // 1..size_of_merge_sets
  // size_of_merge_sets+1...size_of_merge_sets+size_of_new_edges
  std::vector<std::pair<longint, longint> > output;
  output.push_back(std::make_pair(merge_sets_.size(), newEdges_.size()));
  output.insert(output.end(), merge_sets_.begin(), merge_sets_.end());
  output.insert(output.end(), newEdges_.begin(), newEdges_.end());

  mpi::all_gather(*(system_->comm), output, global_merge_sets);

  for (GlobaleMergeSets::iterator gms = global_merge_sets.begin();
       gms != global_merge_sets.end(); ++gms) {
    for (std::vector<std::pair<longint, longint> >::iterator itm = gms->begin();
         itm != gms->end();) {
      longint merge_set_size = itm->first;
      longint new_edge_size = itm->second;
      itm++;
      for (int i = 0; i < merge_set_size; i++, itm++) {
        LOG4ESPP_DEBUG(theLogger, "Merge sets " << itm->first << " with " << itm->second);
        mergeResIdSets(itm->first, itm->second);
      }
      LOG4ESPP_DEBUG(theLogger, "Update new edges. " << new_edge_size);
      // Update new edges.
      for (int i = 0; i < new_edge_size; i++, itm++) {
        newEdge(itm->first, itm->second);
      }
    }
  }
  merge_sets_.clear();
  newEdges_.clear();
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

void TopologyManager::registerPython() {
  using namespace espressopp::python;

  class_<TopologyManager, shared_ptr<TopologyManager>, bases<Extension> >
      ("integrator_TopologyManager", init<shared_ptr<System> >())
      .def("connect", &TopologyManager::connect)
      .def("disconnect", &TopologyManager::disconnect)
      .def("rebuild", &TopologyManager::Rebuild)
      .def("observe_tuple", &TopologyManager::observeTuple)
      .def("observe_triple", &TopologyManager::observeTriple)
      .def("observe_quadruple", &TopologyManager::observeQuadruple)
      .def("initialize", &TopologyManager::InitializeTopology)
      .def("print_topology", &TopologyManager::PrintTopology)
      .def("get_neighbour_lists", &TopologyManager::getNeighbourLists);
}

}  // end namespace integrator
}  // end namespace espressoppp
