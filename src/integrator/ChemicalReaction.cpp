/*
   Copyright (C) 2014-2016
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

#include "ChemicalReaction.hpp"

#include "storage/Storage.hpp"
#include "iterator/CellListIterator.hpp"
#include "esutil/RNG.hpp"
#include "bc/BC.hpp"
#include "storage/NodeGrid.hpp"
#include "storage/DomainDecomposition.hpp"
#include "FixDistances.hpp"

namespace espressopp {
namespace integrator {
LOG4ESPP_LOGGER(ReactionCutoff::theLogger, "ReactionCutoff");
LOG4ESPP_LOGGER(ReactionCutoffStatic::theLogger, "ReactionCutoffStatic");
LOG4ESPP_LOGGER(ReactionCutoffRandom::theLogger, "ReactionCutoffStatic");

void ReactionCutoff::registerPython() {
  using namespace espressopp::python;// NOLINT
  class_<ReactionCutoff, shared_ptr<integrator::ReactionCutoff>, boost::noncopyable>
    ("integrator_ReactionCutoff", no_init);
}

bool ReactionCutoffStatic::check(Particle &p1, Particle &p2) {
  Real3D distance = p1.position() - p2.position();
  real distance_2 = distance.sqr();

  return distance_2 >= min_cutoff_sqr_ && distance_2 < max_cutoff_sqr_;
}

void ReactionCutoffStatic::registerPython() {
  using namespace espressopp::python;// NOLINT

  class_<ReactionCutoffStatic, bases<integrator::ReactionCutoff>,
  boost::shared_ptr<integrator::ReactionCutoffStatic> >
    ("integrator_ReactionCutoffStatic", init<real, real>());
}

bool ReactionCutoffRandom::check(Particle &p1, Particle &p2) {
  real random_cutoff_ = fabs(generator_()) + eq_distance_;
  real random_cutoff_sqr_ = random_cutoff_ * random_cutoff_;
  Real3D distance = p1.position() - p2.position();
  real distance_2 = distance.sqr();

  return distance_2 < random_cutoff_sqr_;
}

void ReactionCutoffRandom::registerPython() {
  using namespace espressopp::python;// NOLINT

  class_<ReactionCutoffRandom, bases<integrator::ReactionCutoff>,
  boost::shared_ptr<integrator::ReactionCutoffRandom> >
    ("integrator_ReactionCutoffRandom", init<real, real, longint>());
}

/** Checks if the particles pair is valid. */
bool Reaction::IsValidPair(Particle &p1, Particle &p2, ParticlePair &particle_order) {
  LOG4ESPP_DEBUG(theLogger, "entering Reaction::IsValidPair");

  if (IsValidState(p1, p2, particle_order)) {
    real W = (*rng_)();
    real p = rate_ * (*dt_) * (*interval_);

    if ((W < p) && reaction_cutoff_->check(p1, p2)) {
      LOG4ESPP_DEBUG(theLogger, "valid pair to bond " << p1.id() << "-" << p2.id());
      return true;
    }
  }

  return false;
}

LOG4ESPP_LOGGER(Reaction::theLogger, "Reaction");

/** Checks if the particles has correct state. */
bool Reaction::IsValidState(Particle &p1, Particle &p2, ParticlePair &correct_order) {
  if ((p1.res_id() == p2.res_id()) && !intramolecular_)
    return false;

  int p1_state = p1.state();
  int p2_state = p2.state();

  // States has to be always positive or zero.
  assert(p1_state >= 0 && p2_state >= 0);

  // Case when both types are the same.
  if ((type_1_ == type_2_) && (p1.type() == type_1_) && (p1.type() == p2.type())) {
    if ((p1_state >= min_state_1_) && (p1_state < max_state_1_) &&
        (p2_state >= min_state_2_) && (p2_state < max_state_2_)) {
      correct_order.first = &p1;
      correct_order.second = &p2;
      return true;
    } else if ((p2_state >= min_state_1_) && (p2_state < max_state_1_) &&
        (p1_state >= min_state_2_) && (p1_state < max_state_2_)) {
      correct_order.first = &p2;
      correct_order.second = &p1;
      return true;
    }
  } else if (type_1_ != type_2_) {// inhomogenious case.
    if (  (p1.type() == type_1_) && (p2.type() == type_2_)
       && ((p1_state >= min_state_1_) && (p1_state < max_state_1_))
       && ((p2_state >= min_state_2_) && (p2_state < max_state_2_))) {
      correct_order.first = &p1;
      correct_order.second = &p2;
      return true;
    } else if (  (p1.type() == type_2_) && (p2.type() == type_1_)
              && ((p1_state >= min_state_2_) && (p1_state < max_state_2_))
              && ((p2_state >= min_state_1_) && (p2_state < max_state_1_))) {
      correct_order.first = &p2;
      correct_order.second = &p1;
      return true;
    }
  }

  return false;
}

bool Reaction::IsValidStateT_1(Particle &p) {
  if (p.type() != type_1_)
    throw std::runtime_error("Particle has wrong type.");

  // States has to be always positive or zero.
  int p_state = p.state();

  assert(p_state >= 0);
  return p_state >= min_state_1_ && p_state < max_state_1_;
}

bool Reaction::IsValidStateT_2(Particle &p) {
  if (p.type() != type_2_)
    throw std::runtime_error("Particle has wrong type.");

  // States has to be always positive or zero.
  int p_state = p.state();

  assert(p_state >= 0);
  return p_state >= min_state_2_ && p_state < max_state_2_;
}

std::set<Particle *> Reaction::PostProcess_T1(Particle &p, Particle &partner) {
  std::set<Particle *> output;
  std::vector<Particle *> ret;

  for (std::vector<shared_ptr<integrator::ChemicalReactionPostProcess> >::iterator it =
      post_process_T1.begin();
      it != post_process_T1.end(); ++it) {
    ret = (*it)->process(p, partner);
    output.insert(ret.begin(), ret.end());
  }

  return output;
}

std::set<Particle *> Reaction::PostProcess_T2(Particle &p, Particle &partner) {
  std::set<Particle *> output;
  std::vector<Particle *> ret;

  for (std::vector<shared_ptr<integrator::ChemicalReactionPostProcess> >::iterator it =
      post_process_T2.begin();
      it != post_process_T2.end(); ++it) {
    ret = (*it)->process(p, partner);
    output.insert(ret.begin(), ret.end());
  }

  return output;
}

void Reaction::registerPython() {
  using namespace espressopp::python;// NOLINT
  class_<Reaction, shared_ptr<integrator::Reaction> >
    ("integrator_Reaction",

      // type_1, type_2, delta_1, delta_2, min_state_1, max_state_1,
      // min_state_2, max_state_2, cutoff, rate, fpl, intramolecular
      init<int, int, int, int, int, int, int, int, real,
      shared_ptr<FixedPairList>, bool>())
    .add_property("type_1", &Reaction::type_1, &Reaction::set_type_1)
    .add_property("type_2", &Reaction::type_2, &Reaction::set_type_2)
    .add_property("delta_1", &Reaction::delta_1, &Reaction::set_delta_1)
    .add_property("min_state_1", &Reaction::min_state_1, &Reaction::set_min_state_1)
    .add_property("max_state_1", &Reaction::max_state_1, &Reaction::set_max_state_1)
    .add_property("delta_2", &Reaction::delta_2, &Reaction::set_delta_2)
    .add_property("min_state_2", &Reaction::min_state_2, &Reaction::set_min_state_2)
    .add_property("max_state_2", &Reaction::max_state_2, &Reaction::set_max_state_2)
    .add_property("rate", &Reaction::rate, &Reaction::set_rate)
    .add_property("intramolecular", &Reaction::intramolecular, &Reaction::set_intramolecular)
    .add_property("active", &Reaction::active, &Reaction::set_active)
    .def("add_postprocess", &Reaction::AddPostProcess)
    .def("set_cutoff", &Reaction::SetReactionCutoff);
}

/** DissociationReaction */
LOG4ESPP_LOGGER(DissociationReaction::theLogger, "DissociationReaction");

/** Checks if the particles pair is valid. */
bool DissociationReaction::IsValidPair(Particle &p1, Particle &p2, ParticlePair &particle_order) {
  LOG4ESPP_DEBUG(theLogger, "entering DissociationReaction::IsValidPair");

  if (IsValidState(p1, p2, particle_order)) {
    real W = (*rng_)();

    if (rate_ > 0.0) {
      Real3D distance;

      bc_->getMinimumImageVectorBox(distance, p1.position(), p2.position());

      real distance_2 = distance.sqr();

      // Break the bond when the distance exceed the cut_off with some probability.
      if ((distance_2 > break_cutoff_sqr_) && (W < rate_ * (*dt_) * (*interval_))) {
        LOG4ESPP_DEBUG(theLogger,
            "Break the bond, " << p1.id() << "-" << p2.id()
                               << " d_2=" << distance_2 << " cutoff_sqr=" << break_cutoff_sqr_);
        return true;
      }
    }

    // Break the bond randomly.
    if (W < diss_rate_ * (*dt_) * (*interval_)) {
      LOG4ESPP_DEBUG(theLogger, "Break the bond randomly " << p1.id() << "-" << p2.id());
      return true;
    }
  }

  return false;
}

void DissociationReaction::registerPython() {
  LOG4ESPP_DEBUG(theLogger, "register dissociation reaction");
  using namespace espressopp::python;// NOLINT
  class_<DissociationReaction, bases<Reaction>, shared_ptr<integrator::DissociationReaction>
  >
    ("integrator_DissociationReaction",

      // type_1, type_2, delta_1, delta_2, min_state_1, max_state_1, min_state_2,
      // max_state_2, break_cutoff, break_rate, fpl
      init<int, int, int, int, int, int, int, int, real, real, shared_ptr
      <FixedPairList> >())
    .add_property("type_1",
      &DissociationReaction::type_1, &DissociationReaction::set_type_1)
    .add_property("type_2",
      &DissociationReaction::type_2, &DissociationReaction::set_type_2)
    .add_property("delta_1",
      &DissociationReaction::delta_1, &DissociationReaction::set_delta_1)
    .add_property("min_state_1",
      &DissociationReaction::min_state_1,
      &DissociationReaction::set_min_state_1)
    .add_property("max_state_1",
      &DissociationReaction::max_state_1,
      &DissociationReaction::set_max_state_1)
    .add_property("delta_2",
      &DissociationReaction::delta_2, &DissociationReaction::set_delta_2)
    .add_property("min_state_2",
      &DissociationReaction::min_state_2,
      &DissociationReaction::set_min_state_2)
    .add_property("intramolecular",
      &DissociationReaction::intramolecular,
      &DissociationReaction::set_intramolecular)
    .add_property("max_state_2",
      &DissociationReaction::max_state_2,
      &DissociationReaction::set_max_state_2)
    .add_property("cutoff", &DissociationReaction::cutoff, &DissociationReaction::set_cutoff)
    .add_property("rate",
      &DissociationReaction::rate, &DissociationReaction::set_rate)
    .add_property("diss_rate",
      &DissociationReaction::diss_rate, &DissociationReaction::set_diss_rate)
    .add_property("active", &DissociationReaction::active, &DissociationReaction::set_active)
    .def("add_postprocess", &DissociationReaction::AddPostProcess);
}
}// namespace integrator
}// namespace espressopp
