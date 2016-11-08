/*
  Copyright (C) 2015
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
#include "ChemicalConversion.hpp"

using namespace espressopp;  //NOLINT

namespace espressopp {
namespace analysis {

real ChemicalConversion::compute_real() const {
  System& system = getSystemRef();
  CellList realCells = system.storage->getRealCells();

  longint local_count = 0;
  for (iterator::CellListIterator cit(realCells); !cit.isDone(); ++cit) {
    if (cit->type() == p_type)
      local_count++;
  }

  longint global_count = 0;
  boost::mpi::all_reduce(*getSystem()->comm, local_count, global_count, std::plus<longint>());

  real value = global_count;
  if (!absolute_value_)
    value = global_count / total_value;
  // Send value via signal.
  onValue(value);

  return value;
}

void ChemicalConversion::registerPython() {
  using namespace espressopp::python;  //NOLINT
  class_<ChemicalConversion, bases<Observable>, boost::noncopyable>
    ("analysis_ChemicalConversion",
        init< shared_ptr<System>, longint, longint >())
        .def(init<shared_ptr<System>, longint>());
}

real ChemicalConversionTypeSequence::compute_real() const {
  System& system = getSystemRef();
  CellList realCells = system.storage->getRealCells();

  longint local_count = 0;
  std::map<longint, longint> local_pid_type;
  for (ParticleGroup::iterator it = particle_group_->begin(); it != particle_group_->end(); ++it) {
    local_pid_type.insert(std::make_pair(it->id(), it->type()));
  }

  // Do the computation on root node as we have to match type sequence and
  // some particles can be splited among nodes.
  std::vector<std::map<longint, longint> > global_pid_type;

  real result = 0.0;
  if (system.comm->rank() == 0) {
    mpi::gather(*(system.comm), local_pid_type, global_pid_type, 0);

    std::map<longint, longint> global_seq;
    for (std::vector<std::map<longint, longint> >::iterator it = global_pid_type.begin();
         it != global_pid_type.end(); ++it ) {
      global_seq.insert(it->begin(), it->end());
    }

    // Flying window, looking for the sequence.
    size_t window_size = type_seq_.size();
    size_t particle_number = global_seq.size();
    if (particle_number % window_size != 0)
      throw std::runtime_error("Number of particles does not corresponds to the sequence length");

    longint counter = 0;

    for (std::map<longint, longint>::iterator it = global_seq.begin(); it != global_seq.end();) {
      bool valid = true;
      for (size_t j = 0; j < window_size; j++, ++it) {
        if (it->second != type_seq_[j] && valid)
          valid = false;
      }
      if (valid)
        counter++;
    }

  } else {
    mpi::gather(*(system.comm), local_pid_type, global_pid_type, 0);
  }

  mpi::broadcast(*(system.comm), result, 0);

  // Send value via signal.
  onValue(result);

  return result;
}

void ChemicalConversionTypeSequence::registerPython() {
  using namespace espressopp::python;  //NOLINT
  class_<ChemicalConversionTypeSequence, bases<Observable>, boost::noncopyable>
      ("analysis_ChemicalConversionTypeSequence",
       init< shared_ptr<System>, shared_ptr<ParticleGroup>, longint >())
      .def("set_sequence", &ChemicalConversionTypeSequence::setSequence);
}


}  // end namespace analysis
}  // end namespace espressopp
