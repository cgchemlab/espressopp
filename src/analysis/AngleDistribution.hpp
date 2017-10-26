/*
  Copyright (C) 2017
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
#ifndef _ANALYSIS_ANGLEDISTRIBUTION_HPP
#define _ANALYSIS_ANGLEDISTRIBUTION_HPP

#include "Observable.hpp"
#include "python.hpp"
#include "types.hpp"
#include "FixedTripleList.hpp"
#include "integrator/TopologyManager.hpp"
#include <set>

namespace espressopp {
namespace analysis {

/** Class to compute the static structure function of the system. */
class AngleDistribution : public Observable {
 public:
  AngleDistribution(shared_ptr<System> system) : Observable(system) {}

  ~AngleDistribution() {}
  virtual real compute() const;

  python::list computeArray(int splitN) const;

  static void registerPython();

 private:
  void registerTriplet(shared_ptr<FixedTripleList> ftl) { ftlList.insert(ftl); }
  void loadFromTopologyManager(shared_ptr<integrator::TopologyManager> topol_man) {
    auto triples = topol_man->getTriples();
    ftlList.insert(triples.begin(), triples.end());
  }

  /**
   * Return angle between p1-p2-p3 particles.
   * @param p1 position of p1 particle.
   * @param p2 position of p2 particle.
   * @param p3 position of p3 particle.
   * @return angle in radians [0, pi]
   */
  real computeAngle(Real3D dist12, Real3D dist32) const;

  std::set<shared_ptr<FixedTripleList> > ftlList;
};

}  // end namespace analysis
}  // end namespace espressopp

#endif
