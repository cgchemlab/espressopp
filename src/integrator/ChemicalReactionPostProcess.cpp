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

/** PostProcessChangeProperty
 *
 * Change property of particles that are involved in the reaction.
 *
 * */
LOG4ESPP_LOGGER(PostProcessChangeProperty::theLogger, "PostProcessChangeProperty");
void PostProcessChangeProperty::addChangeProperty(int type_id,
                                                  boost::shared_ptr<TopologyParticleProperties> new_property) {
  std::pair<TypeParticlePropertiesMap::iterator, bool> ret;


  ret = type_properties_.insert(
      std::pair<int, boost::shared_ptr<TopologyParticleProperties> >(type_id, new_property));

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

  LOG4ESPP_DEBUG(theLogger, "type " << it->second->type());

  if (it != type_properties_.end()) {
    mod = it->second->updateParticleProperties(&p1);
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


/** PostProcessChangePropertyByTopologyManager **/
LOG4ESPP_LOGGER(PostProcessChangePropertyByTopologyManager::theLogger, "PostProcessChangePropertyByTopologyManager");
void PostProcessChangePropertyByTopologyManager::addChangeProperty(
    int type_id, boost::shared_ptr<TopologyParticleProperties> new_property) {
  tm_->registerLocalPropertyChange(type_id, new_property);
}

/** Updates the properties of the particle
 * @param p1 The particle that will be updated
 * @param partner The reference to the partner of p1 particle.
 *
 * */
std::vector<Particle *> PostProcessChangePropertyByTopologyManager::process(Particle &p1, Particle &partner) {
  std::vector<Particle *> mod_particles;
  LOG4ESPP_DEBUG(theLogger, "Entering PostProcessChangeProperty::process()");
  tm_->invokeParticlePropertiesChange(p1.id());
  return mod_particles;
}

void PostProcessChangePropertyByTopologyManager::registerPython() {
  using namespace espressopp::python;// NOLINT

  class_<PostProcessChangePropertyByTopologyManager, bases<integrator::ChemicalReactionPostProcess>,
         boost::shared_ptr<integrator::PostProcessChangePropertyByTopologyManager> >
      ("integrator_PostProcessChangePropertyByTopologyManager", init<shared_ptr<TopologyManager> >())
      .def("add_change_property", &PostProcessChangePropertyByTopologyManager::addChangeProperty);
}

/**
 * This PostProcess will remove bonds from defined bond list (FixedPairList)
 * that are defined between atom involved in pair creation and other atoms.
 * User can defined the number of bonds to remove from this list.
 */
LOG4ESPP_LOGGER(PostProcessRemoveBond::theLogger, "PostProcessRemoveBond");

std::vector<Particle *> PostProcessRemoveBond::process(Particle &p, Particle &partner) {
  LOG4ESPP_DEBUG(theLogger, "Entering PostProcessRemoveBond::operator()");

  std::vector<Particle *> ret;
  // Removes n-bonds.
  if (fpl_->removeByPid1(p.id(), false, false, nr_)) {
    LOG4ESPP_DEBUG(theLogger, "removed bond ");
  }

  return ret;
}

void PostProcessRemoveBond::registerPython() {
  using namespace espressopp::python;// NOLINT

  boost::python::implicitly_convertible<shared_ptr<FixedPairListLambda>, shared_ptr<FixedPairList> >();

  class_<PostProcessRemoveBond, bases<integrator::ChemicalReactionPostProcess>,
         boost::shared_ptr<integrator::PostProcessRemoveBond> >
      ("integrator_PostProcessRemoveBond", init<shared_ptr<FixedPairList>, int>());
}

/** ChangeNeighboursProperty
 *
 * Modyfie property of particles that are neighbour of the new bond pair.
 */
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


/** ChangePropertyOnState
 *
 *  Set the new property of particle whenever the state of particle reaches desired one.
 *  The property are defined in the same way as in PostProcessChangeProperty.
 *
 * */
LOG4ESPP_LOGGER(PostProcessChangePropertyOnState::theLogger, "PostProcessChangePropertyOnState");

std::vector<Particle *> PostProcessChangePropertyOnState::process(Particle &p1, Particle &partner) {
  std::vector<Particle *> ret_val;

  boost::unordered_map<std::pair<longint, longint>, shared_ptr<TopologyParticleProperties> >::const_iterator pp_it =
      type_state_pp_.find(std::make_pair(p1.type(), p1.state()));

  bool mod = false;
  if (pp_it != type_state_pp_.end()) {
    mod = pp_it->second->updateParticleProperties(&p1);
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
      .def("add_change_property", &PostProcessChangePropertyOnState::addChangeProperty);
}


/** PostProcessRemoveNeighbourBond */
void espressopp::integrator::PostProcessRemoveNeighbourBond::registerPython() {
  using namespace espressopp::python;  // NOLINT

  class_<PostProcessRemoveNeighbourBond, bases<ChemicalReactionPostProcess>,
         boost::shared_ptr<PostProcessRemoveNeighbourBond> >
      ("integrator_PostProcessRemoveNeighbourBond", init<shared_ptr<integrator::TopologyManager> >())
      .def("add_bond_to_remove", &PostProcessRemoveNeighbourBond::registerBondToRemove);
}

}// namespace integrator
}// namespace espressopp
