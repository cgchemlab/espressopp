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

#include <vector>
#include "python.hpp"
#include "FixDistances.hpp"

#include "bc/BC.hpp"
#include "types.hpp"
#include "System.hpp"
#include "storage/Storage.hpp"
#include "iterator/CellListIterator.hpp"

namespace espressopp {
using namespace iterator;  //NOLINT
namespace integrator {

LOG4ESPP_LOGGER(FixDistances::theLogger, "FixDistances");

FixDistances::FixDistances(shared_ptr<System> system)
    : Extension(system), has_types_(false) {
  LOG4ESPP_INFO(theLogger, "FixDistances");
  type = Extension::Constraint;
}

FixDistances::FixDistances(shared_ptr<System> system, longint anchor_type, longint target_type)
    : Extension(system), anchor_type_(anchor_type), target_type_(target_type), has_types_(true) {
  LOG4ESPP_INFO(theLogger, "FixDistances");
  type = Extension::Constraint;

  // In theory we would like to observe the type of particles but only if types are provided.
  sigOnParticlesChanged = system->storage->onParticlesChanged.connect
      (boost::bind(&FixDistances::onParticlesChanged, this));
}

void FixDistances::disconnect() {
  aftIntP_.disconnect();
  sigOnParticlesChanged.disconnect();
}

void FixDistances::connect() {
  // Run after position is computed.
  aftIntP_  = integrator->aftIntP.connect(boost::bind(&FixDistances::restore_positions, this));
  if (has_types_) {
    System &system = getSystemRef();
    sigOnParticlesChanged = system.storage->onParticlesChanged.connect
        (boost::bind(&FixDistances::onParticlesChanged, this));
  }
}

void FixDistances::onParticlesChanged() {
  if (!has_types_)
    return;
  std::vector<longint> affected_particles;

  System &system = getSystemRef();
  for (Triplets::iterator it = distance_triplets_.begin(); it != distance_triplets_.end();
       ++it) {
    Particle *anchor = system.storage->lookupLocalParticle(it->first);
    Particle *dst = system.storage->lookupLocalParticle(it->second.first);
    if (anchor && dst) {
      if (anchor->type() != anchor_type_ || dst->type() != target_type_) {
        distance_triplets_.erase(it);
        affected_particles.push_back(it->second.first);
      }
    }
  }

  if (affected_particles.size() > 0 && post_process_) {
    LOG4ESPP_DEBUG(theLogger, "Affected particles " << affected_particles.size());
    for (int i = 0; i < affected_particles.size(); i++) {
      Particle *p1 = system.storage->lookupLocalParticle(affected_particles[i]);
      LOG4ESPP_DEBUG(theLogger, "particle " << p1->id() << " ghost: " << p1->ghost());
      post_process_->operator()(*p1);
    }
  }
}

void FixDistances::restore_positions() {
  System &system = getSystemRef();
  const bc::BC &bc = *(system.bc);
  for (Triplets::iterator it = distance_triplets_.begin(); it != distance_triplets_.end();
       ++it) {
    Particle *anchor = system.storage->lookupLocalParticle(it->first);
    Particle *dst = system.storage->lookupLocalParticle(it->second.first);
    real dist = it->second.second;

    if (anchor != NULL && dst != NULL) {
      Real3D anchor_pos = anchor->position();
      Real3D dst_pos = dst->position();
      Real3D trans;
      bc.getMinimumImageVector(trans, dst_pos, anchor_pos);
      Real3D unit_trans = (1/trans.abs()) * trans;
      Real3D new_trans = dist*unit_trans;

      dst->setPos(anchor_pos + new_trans);
      // Sets the velocity of the anchor and reset the force from anchor.
      dst->setV(anchor->velocity());
      dst->setF((dst->mass()/anchor->mass())*anchor->force());
    }
  }
}

void FixDistances::registerPython() {
  using namespace espressopp::python;  //NOLINT
  class_<FixDistances, shared_ptr<FixDistances>, bases<Extension> >
    ("integrator_FixDistances", init<shared_ptr<System> >())
    .def(init<shared_ptr<System>, int, int>())
    .add_property("size", &FixDistances::size)
    .def("connect", &FixDistances::connect)
    .def("disconnect", &FixDistances::disconnect)
    .def("add_triplet", &FixDistances::add_triplet)
    .def("add_postprocess", &FixDistances::add_postprocess);
}

}  // end namespace integrator
}  // end namespace espressopp
