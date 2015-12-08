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

#include "python.hpp"
#include "DumpTopology.hpp"
#include "storage/Storage.hpp"


namespace espressopp {
namespace io {

void DumpTopology::ObserveTuple(shared_ptr<FixedPairList> fpl) {
  fpl_.push_back(fpl);
  fpl_buffer_.push_back(std::vector<longint>());
}

void DumpTopology::Dump() {
  int fpl_index = 0;
  int current_step = integrator_->getStep();
  for (std::vector<shared_ptr<FixedPairList> >::iterator it = fpl_.begin();
       it != fpl_.end(); ++it, fpl_index++) {
    std::vector<longint> bonds = (*it)->getPairList();
    fpl_buffer_[fpl_index].push_back(current_step);
    fpl_buffer_[fpl_index].push_back(2*bonds.size());
    fpl_buffer_[fpl_index].insert(fpl_buffer_[fpl_index].end(), bonds.begin(), bonds.end());
  }
}

void DumpTopology::ClearBuffer() {
  for (FplBuffer::iterator it = fpl_buffer_.begin(); it != fpl_buffer_.end(); ++it) {
    it->clear();
  }
}

python::list DumpTopology::GetData() {
  python::list ret;
  for (FplBuffer::iterator it = fpl_buffer_.begin(); it != fpl_buffer_.end(); ++it) {
    python::list tpl;
    for (std::vector<longint>::iterator itpl = it->begin(); itpl != it->end(); ++itpl) {
      tpl.append(*itpl);
    }
    ret.append(tpl);
  }
  return ret;
}

// Python wrapping
void DumpTopology::registerPython() {

  using namespace espressopp::python;

  class_<DumpTopology, bases<ParticleAccess>, boost::noncopyable>
      ("io_DumpTopology", init<shared_ptr<System>, shared_ptr<integrator::MDIntegrator> >())
      .def("dump", &DumpTopology::Dump)
      .def("clear_buffer", &DumpTopology::ClearBuffer)
      .def("get_data", &DumpTopology::GetData)
      .def("observe_tuple", &DumpTopology::ObserveTuple);
}

}
}