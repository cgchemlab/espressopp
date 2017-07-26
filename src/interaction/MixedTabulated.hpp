/*
  Copyright (C) 20126
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
#ifndef _INTERACTION_MixedTabulated_HPP
#define _INTERACTION_MixedTabulated_HPP

#include "Potential.hpp"
#include "Interpolation.hpp"
#include "InterpolationLinear.hpp"
#include "InterpolationAkima.hpp"
#include "InterpolationCubic.hpp"
#include "analysis/ChemicalConversion.hpp"

namespace espressopp {
namespace interaction {

class MixedTabulated: public PotentialTemplate<MixedTabulated> {
 public:
  static LOG4ESPP_DECL_LOGGER(theLogger);
  static LOG4ESPP_DECL_LOGGER(theLocalLogger);

  MixedTabulated() {
    setShift(0.0);
    setCutoff(infinity);
    interpolationType = 0;
    mix_value_ = 0.0;
  }

  MixedTabulated(longint itype,
                 const char *filename1,
                 const char *filename2,
                 shared_ptr<analysis::ChemicalConversion> chm,
                 real mix_value,
                 real cutoff) {
    setShift(0.0);
    setCutoff(cutoff);
    interpolationType = itype;
    // Set first table.
    boost::mpi::communicator world;
    switch (itype) {
      case 1:  // Linear interpolation
        table1 = make_shared<InterpolationLinear>();
        table1->read(world, filename1);
        table2 = make_shared<InterpolationLinear>();
        table2->read(world, filename2);
        break;
      case 2:  // Akima interpolation
        table1 = make_shared<InterpolationAkima>();
        table1->read(world, filename1);
        table2 = make_shared<InterpolationAkima>();
        table2->read(world, filename2);
        break;
      case 3:  // Cubic interpolation
        table1 = make_shared<InterpolationCubic>();
        table1->read(world, filename1);
        table2 = make_shared<InterpolationCubic>();
        table2->read(world, filename2);
        break;
      default:
        throw std::runtime_error("wrong itype");
    }
    mix_value_ = mix_value;
    // Connect onValue signal.
    chm->onValue.connect(boost::bind(&MixedTabulated::onValue, this, _1));
  }

  MixedTabulated(longint itype,
                 const char *filename1,
                 const char *filename2,
                 real mix_value,
                 real cutoff) {
    setShift(0.0);
    setCutoff(cutoff);
    interpolationType = itype;
    // Set first table.
    boost::mpi::communicator world;
    switch (itype) {
      case 1:
        table1 = make_shared<InterpolationLinear>();
        table1->read(world, filename1);
        table2 = make_shared<InterpolationLinear>();
        table2->read(world, filename2);
        break;
      case 2:
        table1 = make_shared<InterpolationAkima>();
        table1->read(world, filename1);
        table2 = make_shared<InterpolationAkima>();
        table2->read(world, filename2);
        break;
      case 3:
        table1 = make_shared<InterpolationCubic>();
        table1->read(world, filename1);
        table2 = make_shared<InterpolationCubic>();
        table2->read(world, filename2);
        break;
      default:
        throw std::runtime_error("wrong itype");
    }
    mix_value_ = mix_value;
  }

  /** Returns energy value for given distance square. */
  real _computeEnergySqrRaw(real distSqr) const {
    if (interpolationType != 0) {

      real dist = sqrt(distSqr);
      real e1 = table1->getEnergy(dist);
      real e2 = table2->getEnergy(dist);

      real val = mix_value_ * e1 + (1.0 - mix_value_) * e2;

      LOG4ESPP_DEBUG(theLocalLogger,
                     "Energy, e1=" << e1 << " e2=" << e2
                                   << " x=" << mix_value_ << " val=" << val);
      return val;
    } else {
      return 0.0;
    }
  }

  bool _computeForceRaw(Real3D &force, const Real3D &dist, real distSqr) const {
    if (interpolationType) {
      real ffactor1, ffactor2, ffactor;
      real distrt = sqrt(distSqr);
      ffactor1 = table1->getForce(distrt);
      ffactor2 = table2->getForce(distrt);

      ffactor = mix_value_ * ffactor1 + (1.0 - mix_value_) * ffactor2;

      LOG4ESPP_DEBUG(theLocalLogger,
                     "Force, f1=" << ffactor1 << " f2=" << ffactor2
                                  << " x=" << mix_value_ << " val=" << ffactor);

      force = dist * ffactor;
      return true;
    } else {
      return false;
    }
  }

  static void registerPython();

 private:
  real mix_value() { return mix_value_; }
  void set_mix_value(real s) { mix_value_ = s; }
  longint interpolationType;

  real mix_value_;  //<! The value used for mixing.
  shared_ptr<Interpolation> table1;  //<! Table with U_I
  shared_ptr<Interpolation> table2;  //<! Table with U_II
  void onValue(real value);
};

}  // end namespace interaction
}  // end namespace espressopp
#endif
