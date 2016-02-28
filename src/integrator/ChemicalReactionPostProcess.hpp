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

/** Post process methods to the ChemicalReaction exteions.
 *
 * Each of the espressopp.integrator.Reaction object can accept a list of post-process methods
 * that will perform some actions on the particles after the reaction occures.
 *
 * Reaction: A + B -> A:B
 *
 * The PostProcess method can be defined to act only on A particle, only on B particle or
 * on both particles. Internaly, PostProcess objects are kept on two lists, for A and B particles.
 * First the PostProcess methods are invoked from A list, then from B list. Each time, the process
 * method gets the pointer to A or B particle and pointer to its partner.
 */

/** Base class for PostProcess.**/
class ChemicalReactionPostProcess {
public:
  ChemicalReactionPostProcess() { }

  virtual ~ChemicalReactionPostProcess() { }

  /** The main method that is invoked on the particles that take part in chemical reactions.
   *
   *  We assume that we deal with binary reaction.
   */
  virtual std::vector<Particle *> process(Particle &p, Particle &partner) = 0;

  /** Register this class so it can be used from Python. */
  static void registerPython();

protected:
  shared_ptr<System> system_;
  static LOG4ESPP_DECL_LOGGER(theLogger);
};

/*** PostProcess: change property of particle
 *
 * By using this extension it is possible to modify properties of the particles
 * on-the-fly when reaction occures.
 *
 * Currently it is possible to modify following properties:
 *  - mass
 *  - particle type
 *  - partial charge
 *  - resolution (lambda parameter).
 */
class PostProcessChangeProperty : public integrator::ChemicalReactionPostProcess {
public:
  std::vector<Particle *> process(Particle &p, Particle &partner);

  void AddChangeProperty(int type_id, boost::shared_ptr<ParticleProperties> new_property);
  void RemoveChangeProperty(int type_id);

  /** Register this class so it can be used from Python. */
  static void registerPython();

private:
  TypeParticlePropertiesMap type_properties_;

  static LOG4ESPP_DECL_LOGGER(theLogger);
};

/*** PostProcess: remove bond.
 *
 * By this extension, it is possible to update espressopp.FixedPairList and remove bond(s) from this
 * list. T
 *
 * Example:
 *
 * A(1) + B(2)-B(3) -> A(1)-B(2) + B(3)
 *
 * We have following reaction, molecule A of id 1 forming bond with B of id 2, but B(2) is already
 * bonded with B(3). With this extension, we define which bond list has to be updated, number of
 * bonds that
 * has to be removed.
 */
class PostProcessRemoveBond : public integrator::ChemicalReactionPostProcess {
public:
  PostProcessRemoveBond(shared_ptr<FixedPairList> fpl, int nr):
    fpl_(fpl), nr_(nr) { }

  std::vector<Particle *> process(Particle &p, Particle &partner);

  /** Register this class so it can be used from Python. */
  static void registerPython();

private:
  shared_ptr<FixedPairList> fpl_;

  int nr_;

  static LOG4ESPP_DECL_LOGGER(theLogger);
};
}// namespace integrator
}// namespace espressopp

#endif
