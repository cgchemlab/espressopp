/*
  Copyright (c) 2015
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
#ifndef _IO_DUMPPAIRS_HPP
#define _IO_DUMPPAIRS_HPP

#include "mpi.hpp"
#include <boost/serialization/map.hpp>
#include "types.hpp"
#include "System.hpp"
#include "io/FileBackup.hpp"
#include "ParticleAccess.hpp"
#include "integrator/MDIntegrator.hpp"
#include "storage/Storage.hpp"
#include "iterator/CellListIterator.hpp"
#include "FixedPairList.hpp"

#include "esutil/Error.hpp"

#include <string>

namespace espressopp {
namespace io {
class DumpPairs: public ParticleAccess {
 public:
  DumpPairs(shared_ptr <System> system) : ParticleAccess(system) { }
  ~DumpPairs() {  }

  void perform_action() { dump(); }

  void observeTuple(shared_ptr<FixedPairList>, std::string tuple_name,
                    std::string particle_group = "atoms");
  void dump();

  void update();

  static void registerPython();
 private:
  std::vector<shared_ptr<FixedPairList> > fpl_;
};

}  // end namespace io
}  // end namespace espressopp
#endif
