/*
  Copyright (C) 2018
      Zidan Zhang (zidan.zhang@kuleuven.be)
  
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

#include "python.hpp"
#include "GeneralReactionScheme.hpp"

#include "esutil/RNG.hpp"
#include "bc/BC.hpp"
#include "storage/Storage.hpp"


namespace espressopp {
namespace integrator {

LOG4ESPP_LOGGER(GeneralReactionScheme::theLogger, "GeneralReactionScheme");

GeneralReactionScheme::GeneralReactionScheme(shared_ptr<System> system, int interval)
      : Extension(system), interval_(interval) {
  LOG4ESPP_INFO(theLogger, "GeneralReactionScheme constructed");
  bc_ = getSystem()->bc;
}

void GeneralReactionScheme::disconnect() {
  sig_aftIntV.disconnect();
}

void GeneralReactionScheme::connect() {
  sig_aftIntV = integrator->aftIntV.connect(extensionOrder, boost::bind(&GeneralReactionScheme::updateParticles, this));
}

/**
 * The function get's the particles from CPU and do some operations
 */
void GeneralReactionScheme::updateParticles() {
  if (integrator->getStep() % (interval_) != 0)
    return;

  System &system = getSystemRef();

  // Collect on every CPU the particle properties.
  CellList cl = system.storage->getRealCells();
  longint myN = 0;
  longint totalN = 0;

  std::vector<ParticleData*> particleData;
  particleData.reserve(cl.size());
  for (espressopp::iterator::CellListIterator cit(cl); !cit.isDone(); ++cit) {
    Particle &p = *cit;
    particleData.push_back(new ParticleData(p));
    myN++;
  }

  // Send to master CPU
  longint myRank = system.comm->rank();
  std::vector<std::vector<ParticleData*> > cpus_particleData;
  // Global list of particle data.
  std::vector<ParticleData*> all_particleData;
  all_particleData.reserve(totalN);
  if (myRank == 0) {
    mpi::reduce(*(system.comm), myN, totalN, std::plus<longint>(), 0);
    mpi::gather(*(system.comm), particleData, cpus_particleData, 0);

    for (auto v : cpus_particleData) {
      for (auto p : v) {
        all_particleData.push_back(p);
      }
    }

    // Sort all_particleData by particle id.
    std::sort(all_particleData.begin(), all_particleData.end(), [](ParticleData *p1, ParticleData *p2) {
      return p1->id < p2->id;
    });

    // Get through all pairs.
    for (int i = 0; i < all_particleData.size(); i++) {
      for (int j = i + 1; j < all_particleData.size(); j++) {
        auto p1 = all_particleData[i];
        auto p2 = all_particleData[j];
        Real3D distance;
        bc_->getMinimumImageVectorBox(distance, p1->position, p2->position);
        // absolute distance
        real d = distance.abs();
        // update state + 1
        all_particleData[i]->state++;
        all_particleData[j]->state++;
      }
    }

    // Send back to CPUs
    mpi::broadcast(*(system.comm), all_particleData, 0);
  } else {
    mpi::gather(*(system.comm), particleData, cpus_particleData, 0);
    mpi::broadcast(*(system.comm), all_particleData, 0);
  }

  // Process particles and update locally the data. This happened on all cpus.
  for  (auto particleData : all_particleData) {
    Particle *p = system.storage->lookupLocalParticle(particleData->id);
    if (p) { // if the particle is presented on the cpu, update chemical state and type
      p->setState(particleData->state);
      p->setType(particleData->type);
    }
  }
}

void GeneralReactionScheme::registerPython() {
  using namespace espressopp::python;  // NOLINT

  class_<GeneralReactionScheme, shared_ptr<GeneralReactionScheme>, bases<Extension> >
  ("integrator_GeneralReactionScheme", init<shared_ptr<System>, int>())
      .def("connect", &GeneralReactionScheme::connect)
      .def("disconnect", &GeneralReactionScheme::disconnect);
}

}  // end namespace integrator
}  // end namespace espressopp
