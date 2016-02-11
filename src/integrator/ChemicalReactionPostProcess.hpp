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

// ESPP_CLASS
#ifndef _INTEGRATOR_CHEMICALREACTIONPOSTPROCESS_HPP
#define _INTEGRATOR_CHEMICALREACTIONPOSTPROCESS_HPP

#include <boost/unordered_map.hpp>
#include <boost/signals2.hpp>

#include <utility>
#include <map>
#include <set>
#include <vector>

#include "types.hpp"
#include "logging.hpp"
#include "FixedPairList.hpp"
#include "Particle.hpp"
#include "SystemAccess.hpp"
#include "esutil/ESPPIterator.hpp"

#include "integrator/Extension.hpp"
#include "integrator/VelocityVerlet.hpp"
#include "VerletList.hpp"
#include "interaction/Potential.hpp"


namespace espressopp {
namespace integrator {

typedef std::map<int, boost::shared_ptr<ParticleProperties> > TypeParticlePropertiesMap;

/** PostProcess **/
class ChemicalReactionPostProcess {
 public:
  ChemicalReactionPostProcess() { }
  virtual ~ChemicalReactionPostProcess() { }
  virtual std::vector<Particle*> process(Particle &p, Particle &partner) = 0;

  /** Register this class so it can be used from Python. */
  static void registerPython();

 protected:
  shared_ptr<System> system_;
  static LOG4ESPP_DECL_LOGGER(theLogger);
};


class PostProcessChangeProperty : public integrator::ChemicalReactionPostProcess {
 public:
  std::vector<Particle*> process(Particle &p, Particle &partner);
  void AddChangeProperty(int type_id, boost::shared_ptr<ParticleProperties> new_property);
  void RemoveChangeProperty(int type_id);

  /** Register this class so it can be used from Python. */
  static void registerPython();

 private:
  TypeParticlePropertiesMap type_properties_;
  static LOG4ESPP_DECL_LOGGER(theLogger);
};


class PostProcessRemoveBond : public integrator::ChemicalReactionPostProcess {
 public:
  PostProcessRemoveBond(shared_ptr<FixedPairList> fpl, int nr) : fpl_(fpl), nr_(nr) {}
  std::vector<Particle*> process(Particle &p, Particle &partner);
  void AddPostProcess(shared_ptr<integrator::ChemicalReactionPostProcess> pp) { pp_ = pp; }
  /** Register this class so it can be used from Python. */
  static void registerPython();

 private:
  shared_ptr<FixedPairList> fpl_;
  shared_ptr<integrator::ChemicalReactionPostProcess> pp_;
  int nr_;

  static LOG4ESPP_DECL_LOGGER(theLogger);
};


}  // namespace integrator
}  // namespace espressopp

#endif
