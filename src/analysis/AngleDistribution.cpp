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

#include "AngleDistribution.hpp"
#include <boost/serialization/map.hpp>
#include <cmath>
#include "Configuration.hpp"
#include "FixedTripleList.hpp"
#include "bc/BC.hpp"
#include "esutil/Error.hpp"
#include "iterator/CellListIterator.hpp"
#include "python.hpp"
#include "storage/DomainDecomposition.hpp"

#include <math.h>      // cos and ceil and sqrt
#include <time.h>      // time_t, for particle-distribution-to-cpu time
#include <algorithm>   // std::min
#include <functional>  // std::plus

#ifndef M_PIl
#define M_PIl 3.1415926535897932384626433832795029L
#endif

using namespace espressopp;
using namespace espressopp::iterator;
using namespace std;

namespace espressopp {
namespace analysis {

real AngleDistribution::compute() const { return -1.0; }

python::list AngleDistribution::computeArray(int splitN) const {
  const bc::BC &bc = *getSystemRef().bc;  // boundary conditions
  std::vector<real> angles;
  real histogram[splitN] = {0};
  real dtheta = M_PI / (real)splitN;
  // Iterate over FTL and over angles.
  for (auto ftlIt = ftlList.begin(); ftlIt != ftlList.end(); ftlIt++) {
    for (FixedTripleList::TripleList::Iterator it(**ftlIt); it.isValid(); ++it) {
      Particle &p1 = *it->first;
      Particle &p2 = *it->second;
      Particle &p3 = *it->third;
      Real3D dist12, dist32;
      bc.getMinimumImageVectorBox(dist12, p1.position(), p2.position());
      bc.getMinimumImageVectorBox(dist32, p3.position(), p2.position());
      real angle = computeAngle(dist12, dist32);
      int bin = floor(angle / dtheta);
      histogram[bin] += 1.0;
    }
  }

  real totHistogram[splitN] = {0};
  boost::mpi::all_reduce(*mpiWorld, histogram, splitN, totHistogram, std::plus<real>());

  python::list pyli;
  for (int i = 0; i < splitN; i++) {
    pyli.append(totHistogram[i]);
  }

  return pyli;
}

real AngleDistribution::computeAngle(Real3D dist12, Real3D dist32) const {
  real dist12_sqr = dist12 * dist12;
  real dist32_sqr = dist32 * dist32;
  real dist12_magn = sqrt(dist12_sqr);
  real dist32_magn = sqrt(dist32_sqr);
  real dist1232;
  dist1232 = dist12_magn * dist32_magn;

  real cos_theta = dist12 * dist32 / dist1232;
  if (cos_theta < -1.0)
    cos_theta = -1.0;
  else if (cos_theta > 1.0)
    cos_theta = 1.0;

  return acos(cos_theta);
}

void AngleDistribution::registerPython() {
  using namespace espressopp::python;
  class_<AngleDistribution, bases<Observable> >("analysis_AngleDistribution",
                                                init<shared_ptr<System> >())
      .def("register_triplet", &AngleDistribution::registerTriplet)
      .def("compute", &AngleDistribution::computeArray);
}

}  // end namespace analysis
}  // end namespace espressopp
