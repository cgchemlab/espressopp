/*
  Copyright (C) 2016
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

#include "python.hpp"
#include "FreeOrthorhombicBC.hpp"
#include <cmath>
#include "Real3D.hpp"
#include "Int3D.hpp"
#include "esutil/RNG.hpp"

namespace espressopp {
namespace bc {

LOG4ESPP_LOGGER(FreeOrthorhombicBC::theLogger, "FreeOrthorhombicBC");

/* Constructor */
FreeOrthorhombicBC::FreeOrthorhombicBC(shared_ptr <esutil::RNG> _rng, const Real3D &_boxL, const Int3D &periodicity)
    : BC(_rng) {
  setBoxL(_boxL);
  for (int i = 0; i < 3; i++) {
    if (periodicity[i] == 0)
      periodicity_[i] = false;
    else if (periodicity[i] == 1)
      periodicity_[i] = true;
    else
      std::runtime_error("wrong value");
  }
}

/* Setter method for the box length */
void FreeOrthorhombicBC::setBoxL(const Real3D &_boxL) {
  boxL = _boxL;
  for (int i = 0; i < 3; i++) {
    invBoxL[i] = 1.0 / boxL[i];
    boxL2[i] = 0.5 * boxL[i];
  }
  onBoxDimensionsChanged();
}
void FreeOrthorhombicBC::scaleVolume(real s) {
  boxL *= s;
  boxL2 *= s;
  invBoxL /= s;
  onBoxDimensionsChanged();
}
void FreeOrthorhombicBC::scaleVolume(Real3D s) {
  boxL[0] *= s[0];
  boxL[1] *= s[1];
  boxL[2] *= s[2];
  boxL2[0] *= s[0];
  boxL2[1] *= s[1];
  boxL2[2] *= s[2];
  invBoxL[0] /= s[0];
  invBoxL[1] /= s[1];
  invBoxL[2] /= s[2];
  onBoxDimensionsChanged();
}

/* Returns the minimum image vector between two positions */
void FreeOrthorhombicBC::getMinimumImageVector(Real3D &dist, const Real3D &pos1, const Real3D &pos2) const {
  dist = pos1;
  dist -= pos2;

  for (int i = 0; i < 3; i++) {
    if (periodicity_[i])
      dist[i] -= round(dist[i] * invBoxL[i]) * boxL[i];
  }
}

/* Returns the minimum image vector between two positions */
void FreeOrthorhombicBC::getMinimumImageVectorBox(Real3D &dist, const Real3D &pos1, const Real3D &pos2) const {
  dist = pos1;
  dist -= pos2;

  getMinimumDistance(dist);
}

/* Fold back a nearby position into box */
void FreeOrthorhombicBC::getMinimumDistance(Real3D &dist) const {
  for (int i = 0; i < 3; i++) {
    if (periodicity_[i]) {
      if (dist[i] < -boxL2[i])
        dist[i] += boxL[i];
      else if (dist[i] > boxL2[i])
        dist[i] += boxL[i];
    }
  }
}

/* Returns the minimum image vector between two positions */
void FreeOrthorhombicBC::getMinimumImageVectorX(real dist[3], const real pos1[3], const real pos2[3]) const {
  for (int i = 0; i < 3; i++) {
    dist[i] = pos1[i] - pos2[i];
    if (periodicity_[i])
      dist[i] -= round(dist[i] * invBoxL[i]) * boxL[i];
  }
}

/* Fold an individual coordinate in the specified direction */
void FreeOrthorhombicBC::foldCoordinate(Real3D &pos, Int3D &imageBox, int dir) const {
  if (periodicity_[dir]) {
    int tmp = static_cast<int>(floor(pos[dir] * invBoxL[dir]));

    imageBox[dir] += tmp;
    pos[dir] -= tmp * boxL[dir];

    if (pos[dir] < 0 || pos[dir] >= boxL[dir]) {
      /* slow but safe */
      if (fabs(pos[dir] * invBoxL[dir]) >= INT_MAX / 2) {
        LOG4ESPP_ERROR(theLogger, "folding coordinate: " << pos << " imageBox" << imageBox << " dir: " << dir
            << " result position is out of numerical range. Something wrong!");
        // exit(1);
        // throw std::runtime_error("wrong folded position");
        imageBox[dir] = 0;
        pos[dir] = 0;
      }
    }
  }
}

/* Unfold an individual coordinate in the specified direction */
void FreeOrthorhombicBC::unfoldCoordinate(Real3D &pos, Int3D &imageBox, int dir) const {
  if (periodicity_[dir]) {
    pos[dir] += imageBox[dir] * boxL[dir];
    imageBox[dir] = 0;
  }
}

/* Get random position in the central image box */
void FreeOrthorhombicBC::getRandomPos(Real3D &res) const {
  for (int k = 0; k < 3; k++)
    res[k] = boxL[k];

  res[0] *= (*rng)();
  res[1] *= (*rng)();
  res[2] *= (*rng)();
}

void
FreeOrthorhombicBC::
registerPython() {
  using namespace espressopp::python;
  class_ <FreeOrthorhombicBC, bases<BC>, boost::noncopyable>
  ("bc_FreeOrthorhombicBC", init < shared_ptr < esutil::RNG >, Real3D&, Int3D& > ())
  .add_property("boxL", &FreeOrthorhombicBC::getBoxL, &FreeOrthorhombicBC::setBoxL);
}
}
}
