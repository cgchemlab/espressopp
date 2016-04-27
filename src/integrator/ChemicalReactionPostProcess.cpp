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

#include "ChemicalReactionPostProcess.hpp"

#include "storage/Storage.hpp"
#include "iterator/CellListIterator.hpp"
#include "esutil/RNG.hpp"
#include "bc/BC.hpp"
#include "storage/NodeGrid.hpp"
#include "storage/DomainDecomposition.hpp"

namespace espressopp {
namespace integrator {
LOG4ESPP_LOGGER(ChemicalReactionPostProcess::theLogger, "ChemicalReactionPostProcess");

/** ChemicalReactionPostProcess methods */
void ChemicalReactionPostProcess::registerPython() {
  using namespace espressopp::python;// NOLINT
  class_<ChemicalReactionPostProcess, shared_ptr<integrator::ChemicalReactionPostProcess>,
  boost::noncopyable>
    ("integrator_ChemicalReactionPostProcess", no_init);
}

/** PostProcess - change particle property */
LOG4ESPP_LOGGER(PostProcessChangeProperty::theLogger, "PostProcessChangeProperty");
void PostProcessChangeProperty::addChangeProperty(int type_id,
    boost::shared_ptr<ParticleProperties>             new_property) {
  std::pair<TypeParticlePropertiesMap::iterator, bool> ret;


  ret = type_properties_.insert(
    std::pair<int, boost::shared_ptr<ParticleProperties> >(type_id, new_property));

  if (ret.second == false)
    throw std::runtime_error("Requested type already exists. To replace please remove it firstly");
}

/** Removes change property definition. */
void PostProcessChangeProperty::removeChangeProperty(int type_id) {
  int remove_elements = type_properties_.erase(type_id);


  if (remove_elements == 0) {
    throw std::runtime_error("Invalid type.");
  }
}

/** Updates the properties of the particle
 * @param p1 The particle that will be updated
 * @param partner The reference to the partner of p1 particle.
 *
 * */
std::vector<Particle *> PostProcessChangeProperty::process(Particle &p1, Particle &partner) {
  TypeParticlePropertiesMap::iterator it;
  std::vector<Particle *> mod_particles;


  LOG4ESPP_DEBUG(theLogger, "Entering PostProcessChangeProperty::process()");

  // Process particle p1.
  it = type_properties_.find(p1.type());

  bool mod = false;

  LOG4ESPP_DEBUG(theLogger, "type " << it->second->type);

  if (it != type_properties_.end()) {
    if (it->second->type != NULL) {
      p1.setType(it->second->type);
      mod = true;
    }

    if (it->second->mass != NULL) {
      p1.setMass(it->second->mass);
      mod = true;
    }

    if (it->second->q != NULL) {
      p1.setQ(it->second->q);
      mod = true;
    }

    if (it->second->lambda != NULL) {
      p1.setLambda(it->second->lambda);
      mod = true;
    }

    LOG4ESPP_DEBUG(theLogger, "Modified particle A: " << p1.id());
  }

  if (mod)
    mod_particles.push_back(&(p1));

  return mod_particles;
}

void PostProcessChangeProperty::registerPython() {
  using namespace espressopp::python;// NOLINT

  class_<PostProcessChangeProperty, bases<integrator::ChemicalReactionPostProcess>,
  boost::shared_ptr<integrator::PostProcessChangeProperty> >
    ("integrator_PostProcessChangeProperty", init<>())
    .def("add_change_property", &PostProcessChangeProperty::addChangeProperty)
    .def("remove_change_property", &PostProcessChangeProperty::removeChangeProperty);
}

LOG4ESPP_LOGGER(PostProcessRemoveBond::theLogger, "PostProcessRemoveBond");

std::vector<Particle *> PostProcessRemoveBond::process(Particle &p, Particle &partner) {
  LOG4ESPP_DEBUG(theLogger, "Entering PostProcessRemoveBond::operator()");

  std::vector<Particle *> ret;

  typedef FixedPairList::GlobalPairs GlobalPairs;

  GlobalPairs *gpl = fpl_->getGlobalPairs();
  std::pair<GlobalPairs::const_iterator, GlobalPairs::const_iterator> equalRange = gpl->equal_range(
    partner.id());

  // Removes n-bonds.
  int remove_bonds = 0;

  for (GlobalPairs::const_iterator it = equalRange.first; it != equalRange.second && remove_bonds <
      nr_; ++it) {
    if (fpl_->remove(it->first, it->second)) {
      LOG4ESPP_DEBUG(theLogger, "removed bond "
              << it->first << "-" << it->second << " p=" << p.id() << " partner=" << partner.id());
      remove_bonds++;
    }
  }

  return ret;
}

void PostProcessRemoveBond::registerPython() {
  using namespace espressopp::python;// NOLINT

  class_<PostProcessRemoveBond, bases<integrator::ChemicalReactionPostProcess>,
  boost::shared_ptr<integrator::PostProcessRemoveBond> >
    ("integrator_PostProcessRemoveBond", init<shared_ptr<FixedPairList>, int>());
}


LOG4ESPP_LOGGER(PostProcessChangeNeighboursProperty::theLogger, "PostProcessChangeNeighboursProperty");

std::vector<Particle *> PostProcessChangeNeighboursProperty::process(Particle &p, Particle &partner) {

  topology_manager_->invokeNeighbourPropertyChange(p);

  return std::vector<Particle *>();
}

void PostProcessChangeNeighboursProperty::registerPython() {
  using namespace espressopp::python;// NOLINT

  class_<PostProcessChangeNeighboursProperty, bases<integrator::ChemicalReactionPostProcess>,
         boost::shared_ptr<integrator::PostProcessChangeNeighboursProperty> >
      ("integrator_PostProcessChangeNeighboursProperty", init<shared_ptr<integrator::TopologyManager> >())
      .def("add_change_property", &PostProcessChangeNeighboursProperty::registerNeighbourPropertyChange);

}

LOG4ESPP_LOGGER(PostProcessChangePropertyOnState::theLogger, "PostProcessChangePropertyOnState");

std::vector<Particle *> PostProcessChangePropertyOnState::process(Particle &p1, Particle &partner) {
  std::vector<Particle* > ret_val;

  boost::unordered_map<std::pair<longint, longint>, shared_ptr<ParticleProperties> >::const_iterator pp_it =
      type_state_pp_.find(std::make_pair(p1.type(), p1.state()));

  bool mod = false;
  if (pp_it != type_state_pp_.end()) {
    shared_ptr<ParticleProperties> pp = pp_it->second;
    if (pp->type != NULL) {
      p1.setType(pp->type);
      mod = true;
    }

    if (pp->mass != NULL) {
      p1.setMass(pp->mass);
      mod = true;
    }

    if (pp->q != NULL) {
      p1.setQ(pp->q);
      mod = true;
    }

    if (pp->lambda != NULL) {
      p1.setLambda(pp->lambda);
      mod = true;
    }

    LOG4ESPP_DEBUG(theLogger, "Modified particle id=" << p1.id());
  }

  if (mod)
    ret_val.push_back(&(p1));

  return ret_val;
}

void PostProcessChangePropertyOnState::registerPython() {
  using namespace espressopp::python;  // NOLINT

  class_<PostProcessChangePropertyOnState, bases<integrator::ChemicalReactionPostProcess>,
         boost::shared_ptr<integrator::PostProcessChangePropertyOnState> >
      ("integrator_PostProcessChangePropertyOnState", init<>())
      .def("add_change_property", &PostProcessChangePropertyOnState::addPropertyChange);
}

}// namespace integrator
}// namespace espressopp