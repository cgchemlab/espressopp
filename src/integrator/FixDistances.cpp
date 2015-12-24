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

#include <functional>
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
    : Extension(system) {
  LOG4ESPP_INFO(theLogger, "FixDistances");
  type = Extension::Constraint;
  has_types_ = false;
}

FixDistances::FixDistances(shared_ptr<System> system, longint anchor_type, longint target_type)
    : Extension(system), anchor_type_(anchor_type), target_type_(target_type) {
  LOG4ESPP_INFO(theLogger, "FixDistances");
  type = Extension::Constraint;
  has_types_ = true;
}

void FixDistances::disconnect() {
  aftInitF_.disconnect();
  if (has_types_) {
    aftIntV2_.disconnect();
  }
  sigBeforeSend.disconnect();
  sigAfterRecv.disconnect();
}

void FixDistances::connect() {
  // Calculate force that will move particle, set it after forces are set to zero.
  aftInitF_ = integrator->aftInitF.connect(boost::bind(&FixDistances::restore_positions, this));

  // If use particle types then update constraints at every time step at last.
  System &system = getSystemRef();
  if (has_types_) {
    aftIntV2_ = integrator->aftIntV2.connect(boost::bind(&FixDistances::onParticlesChanged, this));
  }

  sigBeforeSend = system.storage->beforeSendParticles.connect(
      boost::bind(&FixDistances::beforeSendParticles, this, _1, _2));
  sigAfterRecv = system.storage->afterRecvParticles.connect(
      boost::bind(&FixDistances::afterRecvParticles, this, _1, _2));
}

void FixDistances::onParticlesChanged() {
  if (!has_types_)
    return;
  std::vector<Particle*> affected_particles;

  System &system = getSystemRef();
  for (Triplets::iterator it = distance_triplets_.begin(); it != distance_triplets_.end(); ) {
    Particle *anchor = system.storage->lookupRealParticle(it->first);
    Particle *dst = system.storage->lookupRealParticle(it->second.first);
    if (anchor && dst) {
      if (anchor->type() != anchor_type_ || dst->type() != target_type_) {
        affected_particles.push_back(dst);
        it = distance_triplets_.erase(it);
      } else {
        ++it;
      }
    } else {
      ++it;
    }
  }

  // We can do something with the particles that are right now free like changing the type.
  if (affected_particles.size() > 0 && post_process_) {
    LOG4ESPP_DEBUG(theLogger, "Affected particles " << affected_particles.size());
    for (int i = 0; i < affected_particles.size(); i++) {
      Particle *p1 = affected_particles[i];
      LOG4ESPP_DEBUG(theLogger, "particle " << p1->id() << " ghost: " << p1->ghost());
      post_process_->process(*p1);
      // reset force and velocity of released particle.
      p1->setV(Real3D(0.0, 0.0, 0.0));
      p1->setF(Real3D(0.0, 0.0, 0.0));
    }
  }
}

void FixDistances::restore_positions() {
  System &system = getSystemRef();
  const bc::BC &bc = *(system.bc);
  real dt2 = integrator->getTimeStep()*integrator->getTimeStep();
  for (Triplets::iterator it = distance_triplets_.begin(); it != distance_triplets_.end();
       ++it) {
    Particle *anchor = system.storage->lookupRealParticle(it->first);
    Particle *dst = system.storage->lookupLocalParticle(it->second.first);
    real dist = it->second.second;

    if (anchor != NULL && dst != NULL) {
      Real3D anchor_pos = anchor->position();
      Real3D dst_pos = dst->position();
      Real3D unit_trans;
      if (dst_pos.isNaNInf()) {
        LOG4ESPP_ERROR(theLogger, "Particle " << dst->id() << " of anchor " << anchor->id()
            << " has pos: " << dst_pos << " anchor_pos: " << anchor_pos
            << " dst.ghost=" << dst->ghost());
        exit(1);
      }
      // Compute force that will cause to keep dst particle at particular distance from
      // anchor.
      Real3D trans;
      bc.getMinimumImageVector(trans, dst_pos, anchor_pos);
      unit_trans = (1/trans.abs()) * trans;
      Real3D new_trans = (dist - trans.abs())*unit_trans;

      LOG4ESPP_DEBUG(theLogger, "update particle " << dst->id()
          << " new pos: " << anchor_pos + new_trans
          << " unit_trans " << unit_trans
          << " dist " << dist
          << " new_force " << (2*new_trans*(dst->mass())/dt2));

      Real3D newF = (2*new_trans*(dst->mass())/dt2) + dst->getF();
      if (newF.isNaNInf()) {
        std::cout << "new_trans=" << new_trans << std::endl;
        std::cout << "mass=" << dst->mass() << std::endl;
        std::cout << "trans=" << trans.abs() << std::endl;
        std::cout << "dst_pos=" << dst_pos << std::endl;
        std::cout << "anchor_pos=" << anchor_pos << std::endl;
        exit(1);
      }

      dst->setF(newF);
      // dst->setV(anchor->velocity());
      dst->setV(0.0);
    }
  }
}

std::vector<Particle*> FixDistances::release_particle(longint anchor_id, int nr_) {
  int total_size = distance_triplets_.count(anchor_id);
  std::vector<Particle*> mod_particles;
  if (total_size == 0)
    return mod_particles;
  System &system = getSystemRef();

  Particle *p_anchor = system.storage->lookupRealParticle(anchor_id);
  if (!p_anchor) {
    LOG4ESPP_DEBUG(theLogger, "release_particle, anchor_id=" << anchor_id << " not found");
    return mod_particles;
  }

  std::pair<Triplets::iterator, Triplets::iterator> equal_range =
      distance_triplets_.equal_range(anchor_id);
  std::vector<Particle*> tmp;
  int removed = 0;
  for (Triplets::iterator it = equal_range.first; it != equal_range.second; ++it) {
    if (removed == nr_)
      break;
    Particle *p1 = system.storage->lookupLocalParticle(it->second.first);
    if (p1) {
      if (post_process_) {
        tmp = post_process_->process(*p1);
        for (std::vector<Particle*>::iterator it = tmp.begin(); it != tmp.end(); ++it)
          mod_particles.push_back(*it);
      }
      LOG4ESPP_DEBUG(theLogger, "Released particle: " << it->second.first << " of " << anchor_id);
      it = distance_triplets_.erase(it);
      removed++;
      p1->setV(Real3D(0.0, 0.0, 0.0));
      p1->setF(Real3D(0.0, 0.0, 0.0));
    }
  }
  return mod_particles;
}

void FixDistances::add_triplet(longint anchor, longint target, real distance) {
  System &system = getSystemRef();

  // Stores only pairs of real particles.
  Particle *p_anchor = system.storage->lookupRealParticle(anchor);
  Particle *p_target = system.storage->lookupLocalParticle(target);

  bool found = true;
  if (!p_anchor)
    found = false;
  else
    if (!p_target)
      LOG4ESPP_DEBUG(theLogger, "Particle p2 " << target << " not found.");

  if (found) {
    if (has_types_)
      if (has_types_ && (p_anchor->type() != anchor_type_ || p_target->type() != target_type_))
        return;
    distance_triplets_.insert(std::make_pair(anchor, std::pair<longint, real>(target, distance)));
  }
}

void FixDistances::beforeSendParticles(ParticleList &pl, OutBuffer &buf) {
  std::vector<real> toSend;
  int n;
  for (ParticleList::Iterator pit(pl); pit.isValid(); ++pit) {
    longint pid = pit->id();
    n = distance_triplets_.count(pid);
    if (n > 0) {
      std::pair<Triplets::const_iterator, Triplets::const_iterator> equalRange =
        distance_triplets_.equal_range(pid);
      toSend.reserve(toSend.size()+2*n+1);
      toSend.push_back(pid);
      toSend.push_back(2*n);
      for (Triplets::const_iterator it = equalRange.first;
           it != equalRange.second; ++it) {
        toSend.push_back(it->second.first);
        toSend.push_back(it->second.second);
      }
      distance_triplets_.erase(equalRange.first, equalRange.second);
    }
  }
  buf.write(toSend);
}

void FixDistances::afterRecvParticles(ParticleList &pl, InBuffer &buf) {
  std::vector<real> received;
  longint pid1, pid2, n;
  real distance;
  Triplets::iterator it = distance_triplets_.begin();
  buf.read(received);
  int buf_size = received.size();
  int i = 0;
  while (i < buf_size) {
    pid1 = received[i++];
    n = received[i++];
    for (; n > 0; n -= 2) {
      pid2 = received[i++];
      distance = received[i++];
      it = distance_triplets_.insert(
          it, std::make_pair(pid1, std::pair<longint, real>(pid2, distance)));
    }
  }
  if (i != buf_size) {
    LOG4ESPP_ERROR(theLogger,
        "ATTENTION: read garbage during receiving particles\n");
  }
  LOG4ESPP_INFO(theLogger, "received fixed pair list after receive particles");
}

void FixDistances::printTriplets() {
  for (Triplets::iterator it = distance_triplets_.begin(); it != distance_triplets_.end(); ++it) {
     std::cout << it->first << "-" << it->second.first << " d=" <<it->second.second << std::endl;
  }
}

longint FixDistances::totalSize() {
  System &system = getSystemRef();
  longint local_size = distance_triplets_.size();
  longint global_size;
  mpi::all_reduce(*system.comm, local_size, global_size, std::plus<longint>());
  return global_size;
}

void FixDistances::registerPython() {
  using namespace espressopp::python;  //NOLINT
  class_<FixDistances, shared_ptr<FixDistances>, bases<Extension> >
    ("integrator_FixDistances", init<shared_ptr<System> >())
    .def(init<shared_ptr<System>, int, int>())
    .def("get_size", &FixDistances::totalSize)
    .def("connect", &FixDistances::connect)
    .def("disconnect", &FixDistances::disconnect)
    .def("add_triplet", &FixDistances::add_triplet)
    .def("add_postprocess", &FixDistances::add_postprocess)
    .def("print_triplets", &FixDistances::printTriplets);
}

/** Post process after pairs were added.
 *
 * After the reaction is performed, the constraint is removed and dummy particle
 * is released.
 */
LOG4ESPP_LOGGER(PostProcessReleaseParticles::theLogger, "PostProcessReleaseParticles");

std::vector<Particle*> PostProcessReleaseParticles::process(Particle &p) {
  LOG4ESPP_DEBUG(theLogger, "Entering PostProcessReleaseParticles::operator()");
  return fd_->release_particle(p.id(), nr_);
}

void PostProcessReleaseParticles::registerPython() {
  using namespace espressopp::python;  //NOLINT

  class_<PostProcessReleaseParticles, bases<integrator::PostProcess>,
      boost::shared_ptr<integrator::PostProcessReleaseParticles> >
  ("integrator_PostProcessReleaseParticles", init<shared_ptr<integrator::FixDistances>, int>());
}

}  // end namespace integrator
}  // end namespace espressopp
