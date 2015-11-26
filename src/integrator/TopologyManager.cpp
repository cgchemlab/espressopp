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


namespace espressopp {
namespace integrator {

LOG4ESPP_LOGGER(TopologyManager::theLogger, "TopologyManager");

using namespace espressopp::iterator;

TopologyManager::TopologyManager(shared_ptr<System> system) :
    Extension(system), system_(system) {
  LOG4ESPP_INFO(theLogger, "TopologyManager");
  type = Extension::all;
  graph_ = new GraphMap();
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
  aftIntV_ = integrator->aftCalcF.connect(boost::bind(&TopologyManager::exchangeData, this));
}

void TopologyManager::disconnect() {
  aftIntV_.disconnect();
}

void TopologyManager::observeTuple(shared_ptr<FixedPairList> fpl, longint type1, longint type2) {
  fpl->onTupleAdded.connect(
      boost::bind(&TopologyManager::onTupleAdded, this, _1, _2));
  tupleMap_.insert(std::make_pair(std::make_pair(type1, type2), fpl));
}

void TopologyManager::observeTriple(shared_ptr<FixedTripleList> ftl,
                                    longint type1, longint type2, longint type3) {
  tripleMap_.insert(std::make_pair(boost::make_tuple(type1, type2, type3), ftl));
}

void TopologyManager::observeQuadruple(shared_ptr<FixedQuadrupleList> fql, longint type1,
                                       longint type2, longint type3, longint type4) {
  quadrupleMap_.insert(std::make_pair(boost::make_tuple(type1, type2, type3, type4), fql));
}

void TopologyManager::InitializeTopology() {
  // Collect locally the list of edges by iterating over registered tuple lists with bonds.
  typedef std::vector<std::pair<longint, longint> > EdgesVector;
  EdgesVector edges;
  for(TupleMap::iterator it = tupleMap_.begin(); it != tupleMap_.end(); ++it) {
    for (FixedPairList::PairList::Iterator pit(*it->second); pit.isValid(); ++pit) {
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
      newBond(it->first, it->second);
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


void TopologyManager::onTupleAdded(longint pid1, longint pid2) {
  Particle *p1 = system_->storage->lookupRealParticle(pid1);
  Particle *p2 = system_->storage->lookupLocalParticle(pid2);
  if (!p1 || !p2)
    return;
  // Update res_id, if bond joins two molecules with different res_id then merge the sets.
  // p2.res_id() <- p1.res_id()
  if (p1->res_id() != p2->res_id()) {
    merge_sets_.push_back(p1->res_id());
    merge_sets_.push_back(p2->res_id());
  }
  newBond(pid1, pid2);
}

void TopologyManager::newBond(longint pid1, longint pid2) {
  LOG4ESPP_DEBUG(theLogger, "New bond; Adding edge: " << pid1 << "-" << pid2);
  if (graph_->count(pid1) == 0)
    graph_->insert(std::make_pair(pid1, new std::set<longint>()));
  if (graph_->count(pid2) == 0)
    graph_->insert(std::make_pair(pid2, new std::set<longint>()));
  graph_->at(pid1)->insert(pid2);
  graph_->at(pid2)->insert(pid1);

  // Generate angles and dihedrals.
  generateAngles(pid1, pid2);
  generateDihedrals(pid1, pid2);
}

void TopologyManager::generateAngles(longint pid1, longint pid2) {
  std::set<boost::tuple<longint, longint, longint> > triplets;
  std::set<longint> *nb1 = graph_->at(pid1);
  std::set<longint> *nb2 = graph_->at(pid2);
  if (nb1) {
    for (std::set<longint>::iterator it = nb1->begin(); it != nb1->end(); ++it) {
      if (*it == pid1)
        continue;
      triplets.insert(boost::make_tuple(*it, pid1, pid2));
    }
  }
  if (nb2) {
    for (std::set<longint>::iterator it = nb2->begin(); it != nb2->end(); ++it) {
      triplets.insert(boost::make_tuple(pid1, pid2, *it));
    }
  }
}

void TopologyManager::generateDihedrals(longint pid1, longint pid2) {
  std::set<boost::tuple<longint, longint, longint, longint> > quadruplets;
  std::set<longint> *nb1 = graph_->at(pid1);
  std::set<longint> *nb2 = graph_->at(pid2);
  // Case pid2 pid1 <> <>
  if (nb1) {
    for (std::set<longint>::iterator it = nb1->begin(); it != nb1->end(); ++it) {
      if (*it == pid1)
        continue;
      std::set<longint> *nbb1 = graph_->at(*it);
      if (nbb1) {
        for (std::set<longint>::iterator itt = nbb1->begin(); itt != nbb1->end(); ++it) {
          if (*itt == *it)
            continue;
          quadruplets.insert(boost::make_tuple(pid2, pid1, *it, *itt));
        }
      }
    }
  }
  // Case pid1 pid2 <> <>
  if (nb2) {
    for (std::set<longint>::iterator it = nb2->begin(); it != nb2->end(); ++it) {
      if (*it == pid2)
        continue;
      std::set<longint> *nbb2 = graph_->at(*it);
      if (nbb2) {
        for (std::set<longint>::iterator itt = nbb2->begin(); itt != nbb2->end(); ++it) {
          if (*itt == *it)
            continue;
          quadruplets.insert(boost::make_tuple(pid1, pid2, *it, *itt));
        }
      }
    }
  }
  // Case <> pid1 pid2 <>
  if (nb1 && nb2) {
    
  }
}

void TopologyManager::exchangeData() {
  // Collect all message from other CPUs.

  typedef std::vector<std::vector<longint> > GlobaleMergeSets;
  GlobaleMergeSets global_merge_sets;

  mpi::all_gather(*(system_->comm), merge_sets_, global_merge_sets);

  for (GlobaleMergeSets::iterator gms = global_merge_sets.begin();
       gms != global_merge_sets.end(); ++gms) {
    for (std::vector<longint>::iterator itm = gms->begin(); itm != gms->end();) {
      longint set_a_id = *(itm++);
      longint set_b_id = *(itm++);
      mergeResIdSets(set_a_id, set_b_id);
    }
  }
  merge_sets_.clear();
}

void TopologyManager::mergeResIdSets(longint res_id_a, longint res_id_b) {
  PSet *setB = res_particle_ids_[res_id_b];
  // Merge two sets.
  res_particle_ids_[res_id_a]->insert(res_particle_ids_[res_id_b]->begin(),
                                      res_particle_ids_[res_id_b]->end());
  // Update particle res_id;
  for (PSet::iterator it = setB->begin(); it != setB->end(); ++it) {
    Particle *p = system_->storage->lookupLocalParticle(*it);
    if (p) {
      p->setResId(res_id_a);
    }
  }
  res_particle_ids_[res_id_b] = res_particle_ids_[res_id_a];
  delete setB;  // free memory
}

void TopologyManager::Rebuild() {
  for (ResParticleIds::iterator it = res_particle_ids_.begin();
       it != res_particle_ids_.end(); ++it) {
    delete &(*it->second);
  }
  res_particle_ids_.clear();
  // Collect local particle ids.
  CellList cl = system_->storage->getRealCells();
  for (CellListIterator it(cl); it.isValid(); ++it) {
    Particle &p = *it;
    if (res_particle_ids_[p.res_id()] == NULL)
      res_particle_ids_[p.res_id()] = new PSet();
    res_particle_ids_[p.res_id()]->insert(p.id());
  }
  // Sync among CPUs.
  std::vector<ResParticleIds> global_res_particle_ids;
  mpi::all_gather(*(system_->comm), res_particle_ids_, global_res_particle_ids);
  // Update local storage.
  for (std::vector<ResParticleIds>::iterator it = global_res_particle_ids.begin();
       it != global_res_particle_ids.end(); ++it ) {
    for (ResParticleIds::iterator itp = it->begin(); itp != it->end(); ++itp) {
      for (PSet::iterator itps = itp->second->begin(); itps != itp->second->end(); ++itps) {
        if (res_particle_ids_[itp->first] == NULL)
          res_particle_ids_[itp->first] = new PSet();
        res_particle_ids_[itp->first]->insert(*itps);
      }
    }
  }
}

void TopologyManager::registerPython() {
  using namespace espressopp::python;

  class_<TopologyManager, shared_ptr<TopologyManager>, bases<Extension> >
      ("integrator_TopologyManager", init< shared_ptr<System> >())
      .def("connect", &TopologyManager::connect)
      .def("disconnect", &TopologyManager::disconnect)
      .def("rebuild", &TopologyManager::Rebuild)
      .def("observe_tuple", &TopologyManager::observeTuple)
      .def("observe_triple", &TopologyManager::observeTriple)
      .def("observe_quadruple", &TopologyManager::observeQuadruple)
      .def("initialize", &TopologyManager::InitializeTopology)
      .def("print_topology", &TopologyManager::PrintTopology)
      ;
}

}  // end namespace integrator
}  // end namespace espressoppp