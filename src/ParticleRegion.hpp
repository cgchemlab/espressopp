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
#include "integrator/MDIntegrator.hpp"

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
  ParticleRegion(shared_ptr<storage::Storage> _storage, shared_ptr<integrator::MDIntegrator>);
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

  void set_v(real vx, real vy, real vz, bool left_right) {
    if (left_right)
      velocity_left_ = Real3D(vx, vy, vz);
    else
      velocity_right_ = Real3D(vx, vy, vz);
  }

  python::tuple get_v(bool left_right) {
    if (left_right)
      return python::make_tuple(velocity_left_[0], velocity_left_[1], velocity_left_[2]);
    else
      return python::make_tuple(velocity_right_[0], velocity_right_[1], velocity_right_[2]);
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

  // Velocity of modyfication
  Real3D velocity_left_;
  Real3D velocity_right_;

  // pointer to storage object
  shared_ptr<storage::Storage> storage;

  // pointer to integrator
  shared_ptr<integrator::MDIntegrator> integrator_;

  // some signalling stuff to keep track of the particles in cell
  boost::signals2::connection con_changed, sig_aftIntV1, sig_aftIntV2;

  void onParticlesChanged();

  static LOG4ESPP_DECL_LOGGER(theLogger);

 private:
  Real3D left_bottom_;
  Real3D right_top_;

  void updateRegion();
  python::tuple getRegion();
};

}

#endif	/* ParticleRegion_H */

