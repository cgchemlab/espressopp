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

// ESPP_CLASS
#ifndef   __ESPRESSOPP_PARTICLEREGION_H
#define   __ESPRESSOPP_PARTICLEREGION_H

#include "Particle.hpp"
#include "ParticleGroup.hpp"
#include "log4espp.hpp"
#include "types.hpp"
#include <set>
#include <boost/signals2.hpp>
#include "iterator/CellListIterator.hpp"

namespace espressopp {

/**
 * \brief Group of articles in given geometrical region.
 *
 * This is the list of particles that are in given geometrical region.
 * Currently only block shape is supported.
 *
 * The list of particles is updated whenever the particle leave or enter region.
 * Moreover, the region can contains only particles of given type.
 *
 * The ParticleRegion behaves in the same way as ParticleGroup and shares the same interface.
 *
 */
class ParticleRegion : public ParticleGroup {
 public:
  ParticleRegion(shared_ptr <storage::Storage> _storage);
  ~ParticleRegion();

  /**
   * Define the block region by size numbers.
   * @param lo position of left bottom front vertex
   * @param hi position of right top rear vertex.
   */
  void defineRegion(Real3D lo, Real3D hi) {
    left_bottom_ = lo;
    right_top_ = hi;
  }

  void addTypeId(longint type_id) {
    types_.insert(type_id);
    has_types_ = true;
  }

  void removeTypeId(longint type_id) {
    types_.erase(type_id);
    has_types_ = types_.size() > 0;
  }

  // for debugging purpose
  void print();

  bool has(longint pid);

  longint size() { return active.size(); }

  /**
  * \brief begin iterator for active particles
  *
  * @return begin iterator
  */
  iterator begin() {
    return active.begin();
  }

  /**
   * \brief end iterator for active particles
   *
   * @return end iterator
   */
  iterator end() {
    return active.end();
  }

  python::list getParticleIDs();

  static void registerPython();

 protected:
  // list of active particles
  // key: particle id, value: reference to local particle
  std::map<longint, Particle *> active;

  // keep the list of particle types, if empty then all particle types are used.
  std::set<longint> types_;
  bool has_types_;

  // pointer to storage object
  shared_ptr <storage::Storage> storage;

  // some signalling stuff to keep track of the particles in cell
  boost::signals2::connection con_send, con_recv, con_changed;

  void beforeSendParticles(ParticleList &pl, class OutBuffer &buf);
  void afterRecvParticles(ParticleList &pl, class InBuffer &buf);
  void onParticlesChanged();

  static LOG4ESPP_DECL_LOGGER(theLogger);

 private:
  Real3D left_bottom_;
  Real3D right_top_;
};

}

#endif	/* ParticleRegion_H */

