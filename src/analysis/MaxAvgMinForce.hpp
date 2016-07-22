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

// ESPP_CLASS
#ifndef _ANALYSIS_MAXAVGMINFORCE_HPP
#define _ANALYSIS_MAXAVGMINFORCE_HPP

#include "types.hpp"
#include "Observable.hpp"

namespace espressopp {
namespace analysis {

/** Class to get the number of particles in the system. */
class MaxForce : public Observable {
 public:
  MaxForce(shared_ptr<System> system) : Observable(system) {}
  virtual ~MaxForce() {}
  virtual real compute() const;
 private:
  static void registerPython();
};

}
}
#endif
