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
}

TopologyManager::~TopologyManager() {
  disconnect();
}

void TopologyManager::connect() {
  aftIntV_ = integrator->aftCalcF.connect(boost::bind(&TopologyManager::exchangeData, this));
}

void TopologyManager::disconnect() {
  aftIntV_.disconnect();
}

void TopologyManager::observeTuple(shared_ptr<FixedPairList> fpl) {
  fpl->onTupleAdded.connect(
      boost::bind(&TopologyManager::onTupleAdded, this, _1, _2));
  // FIXME: what to do when tuple will be removed?
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
  for (ResParticleIds::iterator it = res_particle_ids_.begin(); it != res_particle_ids_.end(); ++it) {
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
      .def("observe", &TopologyManager::observeTuple)
      ;
}



}  // end namespace integrator
}  // end namespace espressoppp