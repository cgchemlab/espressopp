/*
  Copyright (C)
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

// ESPP_CLASS
#ifndef _INTEGRATOR_FIXDISTANCE_HPP
#define _INTEGRATOR_FIXDISTANCE_HPP

#include <utility>
#include <vector>
#include "types.hpp"
#include "logging.hpp"
#include "Extension.hpp"
#include "ParticleGroup.hpp"
#include "boost/signals2.hpp"
#include "boost/unordered_map.hpp"
#include "Real3D.hpp"
#include "Particle.hpp"
#include "ChemicalReaction.hpp"

namespace espressopp {
namespace integrator {

class FixDistances : public Extension {
 public:
  typedef boost::unordered_multimap<longint, std::pair<longint, real> > Triplets;

  explicit FixDistances(shared_ptr<System> _system);
  FixDistances(shared_ptr<System> _system, longint anchor_type, longint target_type);

  void add_triplet(longint anchor, longint target, real distance) {
    distance_triplets_.insert(std::make_pair(anchor, std::pair<longint, real>(target, distance)));
  }

  void add_postprocess(const shared_ptr<integrator::PostProcessChangesProperty> pp) {
    post_process_ = pp;
  }

  void restore_positions();

  /** Register this class so it can be used from Python. */
  static void registerPython();

 private:
  boost::signals2::connection befIntP_, aftIntP_, sigOnParticlesChanged;
  Triplets distance_triplets_;

  longint anchor_type_, target_type_;
  bool has_types_;

  void connect();
  void disconnect();

  void onParticlesChanged();

  shared_ptr<integrator::PostProcessChangesProperty> post_process_;

  /** Logger */
  static LOG4ESPP_DECL_LOGGER(theLogger);
};

}  // end namespace integrator
}  // end namespace espressopp
#endif
