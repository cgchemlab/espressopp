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

// ESPP_CLASS
#ifndef _INTEGRATOR_FIXDISTANCE_HPP
#define _INTEGRATOR_FIXDISTANCE_HPP

#include <utility>
#include <vector>
#include "types.hpp"
#include "logging.hpp"
#include "integrator/Extension.hpp"
#include "ParticleGroup.hpp"
#include "boost/signals2.hpp"
#include "boost/unordered_map.hpp"
#include "Real3D.hpp"
#include "Particle.hpp"
#include "ChemicalReactionPostProcess.hpp"
#include "analysis/Observable.hpp"

namespace espressopp {
namespace integrator {

/**
 * Distance constraint on particles. The target particle will always be kept on certain
 * distance from anchor particle. Optional it can observer also types of particles, whenever
 * anchor particle or target particle changes the type, the constraint will be released.
 *
 * @param System: The espressopp.System object.
 * @param anchor_type: The type of anchor particle.
 * @param target_type: The type of target particle.
 */
class FixDistances : public Extension {
 public:
  typedef boost::unordered_multimap<longint, std::pair<longint, real> > Triplets;

  explicit FixDistances(shared_ptr<System> _system);
  FixDistances(shared_ptr<System> _system, longint anchor_type, longint target_type);

  /**
   *  Adds triplets of anchor, target and distance.
   *
   *  @param anchor: The id of anchor particle.
   *  @param target: The id of target particle.
   *  @param distance: The distance on which particles will be kept.
   */
  void add_triplet(longint anchor, longint target, real distance, bool force);

  void add_postprocess(const shared_ptr<integrator::PostProcessChangeProperty> pp) {
    post_process_ = pp;
  }

  int totalSize() const;
  void restore_positions();
  std::vector<Particle*> release_particle(longint anchor_id, int nr_=1);

  /** Register this class so it can be used from Python. */
  static void registerPython();

 private:
  boost::signals2::connection aftInitF_, sigOnParticlesChanged, aftIntV_;
  boost::signals2::connection sigBeforeSend, sigAfterRecv;
  Triplets distance_triplets_;

  longint anchor_type_, target_type_;
  bool has_types_;

  shared_ptr<System> system_;

  python::list getAllTriplets();

  void connect();
  void disconnect();

  void onAftIntV();
  void beforeSendParticles(ParticleList& pl, OutBuffer& buf);
  void afterRecvParticles(ParticleList &pl, InBuffer& buf);
  shared_ptr<integrator::PostProcessChangeProperty> post_process_;
  void printTriplets();
  /** Logger */
  static LOG4ESPP_DECL_LOGGER(theLogger);

  /** Timers */
  real timeRestorePosition;
  real timeUpdateList;
  real timeComm;

  esutil::WallTimer wallTimer;
  python::list getTimers();
};

/**
 * PostProcess action. Invoked when target particle is realesd from the constraint.
 *
 * @params fd espressopp.integrator.FixDistances object.
 * @params nr The integer with the number of constraints to release at once.
 */
class PostProcessReleaseParticles : public integrator::ChemicalReactionPostProcess {
 public:
  PostProcessReleaseParticles(shared_ptr<FixDistances> fd, int nr) : fd_(fd), nr_(nr) {}
  std::vector<Particle*> process(Particle &p, Particle &partner);

  longint getOrder() { return 5; }
  std::string getName() { return "PostProcessReleaseParticles"; }

  static void registerPython();
 private:
  int nr_;
  shared_ptr<integrator::FixDistances> fd_;
  /** Logger */
  static LOG4ESPP_DECL_LOGGER(theLogger);
};

class PostProcessJoinParticles : public integrator::ChemicalReactionPostProcess {
 public:
  PostProcessJoinParticles(shared_ptr<FixDistances> fd, real distance) : fd_(fd), distance_(distance) {}
  std::vector<Particle*> process(Particle &p, Particle &partner);

  longint getOrder() { return 5; }
  std::string getName() { return "PostProcessJoinParticles"; }

  static void registerPython();

 private:
  shared_ptr<integrator::FixDistances> fd_;
  real distance_;

  static LOG4ESPP_DECL_LOGGER(theLogger);
};

}  // end namespace integrator
}  // end namespace espressopp
#endif
