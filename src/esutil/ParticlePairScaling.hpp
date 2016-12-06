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

#ifndef _ESUTIL_PARTICLEPAIRSCALING_HPP
#define _ESUTIL_PARTICLEPAIRSCALING_HPP

#include "types.hpp"
#include "VerletList.hpp"
#include "boost/unordered_map.hpp"

namespace espressopp {
namespace esutil {

class ParticlePairScaling {
 public:
  ParticlePairScaling() {
    default_scale_ = 1.0;
    d_scale_factor_ = 0.0;
  }

  ParticlePairScaling(real default_scale, real d_scale_factor, shared_ptr<VerletList> vl, shared_ptr<integrator::MDIntegrator> integrator)
      : default_scale_(default_scale), d_scale_factor_(d_scale_factor) {
    sigOnPairUnexclude = vl->onPairUnexclude.connect(boost::bind(&ParticlePairScaling::addParticlePair, this, _1, _2));
    sigOnPairExclude = vl->onPairExclude.connect(boost::bind(&ParticlePairScaling::deleteParticlePair, this, _1, _2));

    sigOnAftIntV = integrator->aftIntV.connect(boost::bind(&ParticlePairScaling::incrementAllScaleFactors, this));
  }

  ~ParticlePairScaling() {
    sigOnPairUnexclude.disconnect();
    sigOnPairExclude.disconnect();
    sigOnAftIntV.disconnect();
  }

  real getPairScaling(longint pid1, longint pid2);

  static void registerPython();
 private:
  typedef boost::unordered_map<std::pair<longint, longint>, real> ParticlePairScale;
  void addParticlePair(longint pid1, longint pid2);
  void deleteParticlePair(longint pid1, longint pid2);

  void incrementAllScaleFactors();

  ParticlePairScale particle_pair_scale_;
  real default_scale_;
  real d_scale_factor_;
  boost::signals2::connection sigOnPairUnexclude, sigOnPairExclude, sigOnAftIntV;

  static LOG4ESPP_DECL_LOGGER(theLogger);
};

}  // end namespace esutil
}  // end namespace espressopp
#endif
