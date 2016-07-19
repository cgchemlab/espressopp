/*
  Copyright (C) 2016
      Jakub Krajniak (c) (jkrajniak at gmail.com)
  
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
#include "ParticleRegion.hpp"
#include "storage/Storage.hpp"
#include "iterator/CellListIterator.hpp"

namespace espressopp {

LOG4ESPP_LOGGER(ParticleRegion::theLogger, "ParticleRegion");

ParticleRegion::ParticleRegion(shared_ptr<storage::Storage> _storage) : storage(_storage), ParticleGroup() {

  con_changed = storage->onParticlesChanged.connect
      (boost::bind(&ParticleRegion::onParticlesChanged, this));
  has_types_ = false;
}

ParticleRegion::~ParticleRegion() {
  con_changed.disconnect();
}

bool ParticleRegion::has(longint pid) {
  return particles.find(pid) != particles.end();
}

// for debugging purpose
void ParticleRegion::print() {
  std::cout << "####### I have " << active.size() << " active particles";
  std::cout << " region: " << left_bottom_ << " to " << right_top_ << std::endl;
  for(iterator i=begin(); i!=end(); ++i ) {
    std::cout << "pid: " << i->id() << " " << i->type() << std::endl;
  }
}

python::list ParticleRegion::getParticleIDs() {
  python::list particle_ids;
  for(iterator i=begin(); i!=end(); ++i ) {
    particle_ids.append(i->id());
  }

  return particle_ids;
}

void ParticleRegion::onParticlesChanged() {
  LOG4ESPP_DEBUG(theLogger, "ParticleRegion::onParticlesChanges");
  LOG4ESPP_DEBUG(theLogger, "left: " << left_bottom_ << " right: " << right_top_);
  active.clear();

  // Update active list.
  CellList cl = storage->getLocalCells();
  for(espressopp::iterator::CellListIterator cit(cl); !cit.isDone(); ++cit) {
    Particle &p = *cit;
    LOG4ESPP_DEBUG(theLogger, "particle id=" << p.id() << " t= " << p.type() << " p= " << p.position() << " g=" << p.ghost());
    if (!has_types_ || (has_types_ && types_.count(p.type()) == 1)) {
      Real3D pos = p.position();
      bool in_region = true;
      for (int d = 0; d < 3; d++) {
        in_region = in_region && (pos[d] > left_bottom_[d]) && (pos[d] < right_top_[d]);
      }
      if (in_region) {
        LOG4ESPP_DEBUG(theLogger, "insert  p " << p.id() << " in=" << in_region);
        active.insert(std::make_pair(p.id(), &p));
      }
    }
  }
}

void ParticleRegion::registerPython() {
  using namespace espressopp::python;

  class_<ParticleRegion, shared_ptr<ParticleRegion>, bases<ParticleGroup> >
      ("ParticleRegion", init<shared_ptr<storage::Storage> >())
        .def("show", &ParticleRegion::print)
        .def("has", &ParticleRegion::has)
        .def("define_region", &ParticleRegion::defineRegion)
        .def("add_type_id", &ParticleRegion::addTypeId)
        .def("remove_type_id", &ParticleRegion::removeTypeId)
        .def("get_particle_ids", &ParticleRegion::getParticleIDs)
        .def("size", &ParticleRegion::size);
}

}  // end namespace espressopp