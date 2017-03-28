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
#include "ParticlePairScaling.hpp"

namespace espressopp {
namespace esutil {

LOG4ESPP_LOGGER(ParticlePairScaling::theLogger, "ParticlePairScaling");

void ParticlePairScaling::addParticlePair(longint pid1, longint pid2) {
  LOG4ESPP_DEBUG(theLogger, "new pair " << pid1 << "-" << pid2 << " scale: " << default_scale_);
  particle_pair_scale_.insert(std::make_pair(std::make_pair(pid1, pid2), default_scale_));
  particle_pair_scale_.insert(std::make_pair(std::make_pair(pid2, pid1), default_scale_));
}

void ParticlePairScaling::deleteParticlePair(longint pid1, longint pid2) {
  LOG4ESPP_DEBUG(theLogger, "delete pair " << pid1 << "-" << pid2);
  particle_pair_scale_.erase(std::make_pair(pid1, pid2));
  particle_pair_scale_.erase(std::make_pair(pid2, pid1));
}

void ParticlePairScaling::registerPython() {
  using namespace espressopp::python;

  class_<ParticlePairScaling, shared_ptr<ParticlePairScaling>, boost::noncopyable>
      ("esutil_ParticlePairScaling", init<real, real, shared_ptr<VerletList>, shared_ptr<integrator::MDIntegrator> >());
}

real ParticlePairScaling::getPairScaling(longint pid1, longint pid2) {
  ParticlePairScale::const_iterator it = particle_pair_scale_.find(std::make_pair(pid1, pid2));
  if (it == particle_pair_scale_.end()) {
    return 1.0;
  } else {
    return it->second;
  }
}

void ParticlePairScaling::incrementAllScaleFactors() {
  for (ParticlePairScale::iterator it = particle_pair_scale_.begin(); it != particle_pair_scale_.end();) {
    it->second += d_scale_factor_;
    if (it->second >= 1.0) {
      LOG4ESPP_DEBUG(theLogger,
        "scale factor for pair " << it->first.first << "-" << it->first.second << " reached 1.0, pair removed");
      it = particle_pair_scale_.erase(it);
    } else {
      ++it;
    }
  }
}

}
}