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
#include "System.hpp"
#include "ParticleAccess.hpp"
#include "integrator/MDIntegrator.hpp"
#include "storage/Storage.hpp"
#include "FixedPairList.hpp"
#include <vector>

#include "esutil/Error.hpp"

#include <string>

namespace espressopp {
namespace io {
class DumpTopology: public ParticleAccess {
 public:
  DumpTopology(shared_ptr<System> system, shared_ptr<integrator::MDIntegrator> integrator)
      : ParticleAccess(system), integrator_(integrator) { }
  ~DumpTopology() {  }

  void perform_action() { Dump(); }

  void ObserveTuple(shared_ptr<FixedPairList> fpl);
  void Dump();

  python::list GetData();

  static void registerPython();
 private:
  void ClearBuffer();

  shared_ptr<integrator::MDIntegrator> integrator_;
  std::vector<shared_ptr<FixedPairList> > fpl_;
  // Buffer for fpl data.
  typedef std::vector<std::vector<longint> > FplBuffer;
  FplBuffer fpl_buffer_;
};

}  // end namespace io
}  // end namespace espressopp
#endif
