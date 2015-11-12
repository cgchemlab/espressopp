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


namespace espressopp {
namespace integrator {

using namespace espressopp::iterator;

TopologyManager::TopologyManager(shared_ptr<System> system) :
    Extension(system), system_(system) {
  LOG4ESPP_INFO(theLogger, "TopologyManager");
  type = Extension::all;
  is_dirty_ = false;
}

void TopologyManager::observeTuple(shared_ptr<FixedPairList> fpl) {
  fpl->onTupleAdded.connect(
      boost::bind(&TopologyManager::onTupleAdded, this, _1, _2));
}

void TopologyManager::onTupleAdded(longint pid1, longint pid2) {
  Particle *p1 = system_->storage->lookupLocalParticle(pid1);
  Particle *p2 = system_->storage->lookupLocalParticle(pid2);
  if (!p1 || !p2)
    return;
  // Update res_id, if bond joins two molecules with different res_id then merge them.
  if (p1->type() != p2->type()) {

  }
}

void TopologyManager::rebuild() {
  res_particle_ids_.clear();
  CellList cl = system_->storage->getRealCells();
  for (CellListIterator it(cl); it.isValid(); ++it) {
    Particle &p = *it;
    res_particle_ids_.insert(std::make_pair(p.res_id(), p.id()));
  }
  // Sync among CPUs.
  std::vector<ResParticleIds> global_res_particle_ids;
  mpi::all_gather(*(system_->comm), res_particle_ids_, global_res_particle_ids);
}

}  // end namespace integrator
}  // end namespace espressoppp