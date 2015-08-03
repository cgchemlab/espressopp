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

#include <utility>
#include <vector>
#include "python.hpp"
#include "FixDistances.hpp"

#include "bc/BC.hpp"
#include "esutil/RNG.hpp"
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

  // We can do something with the particles that are right now free like chaning the type.
  if (affected_particles.size() > 0 && post_process_) {
    LOG4ESPP_DEBUG(theLogger, "Affected particles " << affected_particles.size());
    for (int i = 0; i < affected_particles.size(); i++) {
      Particle *p1 = system.storage->lookupLocalParticle(affected_particles[i]);
      LOG4ESPP_DEBUG(theLogger, "particle " << p1->id() << " ghost: " << p1->ghost());
      post_process_->process(*p1);
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
      Real3D unit_trans;
      if (dst_pos.isNaNInf()) {
        // Problem with position. Reset it by picking random point on sphere
        LOG4ESPP_DEBUG(theLogger, "Particle " << dst->id() << " of anchor " << anchor->id()
            << " has pos: " << dst_pos << " anchor_pos: " << anchor_pos);
        std::cout << "Problem with position of point " << dst->id() << std::endl;
        unit_trans = system.rng->uniformOnSphere();
      } else {
        Real3D trans;
        bc.getMinimumImageVector(trans, dst_pos, anchor_pos);
        unit_trans = (1/trans.abs()) * trans;
      }
      Real3D new_trans = dist*unit_trans;
      LOG4ESPP_DEBUG(theLogger, "update particle " << dst->id()
          << " new pos: " << anchor_pos + new_trans
          << " unit_trans " << unit_trans
          << " dist " << dist
          << " new_force " << (dst->mass()/anchor->mass())*anchor->force());
      // Sets position of the anchor + new translation and reset the force
      // to the same as is acting on the anchor with scaled valued.
      dst->setPos(anchor_pos + new_trans);
      dst->setF((dst->mass()/anchor->mass())*anchor->force());
      dst->setV(Real3D(0.0, 0.0, 0.0));
    }
  }
}

std::vector<Particle*> FixDistances::release_particle(longint anchor_id) {
  Triplets::iterator fparticle = distance_triplets_.find(anchor_id);

  std::vector<Particle*> mod_particles;
  if (fparticle != distance_triplets_.end()) {
    System &system = getSystemRef();
    Particle *p1 = system.storage->lookupLocalParticle(fparticle->second.first);
    if (p1 != NULL) {
      std::cout << "release particle " << p1->id() << std::endl;
    }
    if (post_process_ && post_process_->process(*p1))
      mod_particles.push_back(p1);
    p1->setV(Real3D(0.0, 0.0, 0.0));
    p1->setF(Real3D(0.0, 0.0, 0.0));
  }

  return mod_particles;

  /*
  std::vector<Particle*> released_particles;

  for (Triplets::iterator it = equal_range.first; it != equal_range.second; ++it) {
    Particle *p1 = system.storage->lookupLocalParticle(it->second.first);
    if (p1 != NULL) {
      std::cout << "release particle " << p1->id() << std::endl;
      released_particles.push_back(p1);
      break;  // release only one at the time;
    }
  }
  distance_triplets_.erase(equal_range.first, equal_range.second);

  if (released_particles.size() > 0) {
    LOG4ESPP_DEBUG(theLogger, "Affected particles " << released_particles.size());
    for (int i = 0; i < released_particles.size(); i++) {
      Particle *p1 = released_particles[i];
      if (p1 != NULL) {
        LOG4ESPP_DEBUG(theLogger, "particle " << p1->id() << "ghost: " << p1->ghost());
        if (post_process_ && post_process_->process(*p1))
          mod_particles.push_back(p1);

        // Reset velocity and force on copartner
        p1->setV(Real3D(0.0, 0.0, 0.0));
        p1->setF(Real3D(0.0, 0.0, 0.0));
      }
    }
  }
  return mod_particles;*/
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

/** Post process after pairs were added.
 *
 * After the reaction is performed, the constraint is removed and dummy particle
 * is released.
 */
LOG4ESPP_LOGGER(PostProcessReleaseParticles::theLogger, "PostProcessReleaseParticles");

std::vector<Particle*> PostProcessReleaseParticles::process(Particle &p1, Particle &p2) {
  LOG4ESPP_DEBUG(theLogger, "Entering PostProcessReleaseParticles::operator()");
  if (nr_ == 1)
    return fd_->release_particle(p1.id());
  return std::vector<Particle*>();
}

void PostProcessReleaseParticles::registerPython() {
  using namespace espressopp::python;  //NOLINT

  class_<PostProcessReleaseParticles, bases<integrator::PostProcess>,
      boost::shared_ptr<integrator::PostProcessReleaseParticles> >
  ("integrator_PostProcessReleaseParticles", init<shared_ptr<integrator::FixDistances>, int>());
}

}  // end namespace integrator
}  // end namespace espressopp
