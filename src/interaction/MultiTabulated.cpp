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
#include "MultiTabulated.hpp"
#include "InterpolationLinear.hpp"
#include "InterpolationAkima.hpp"
#include "InterpolationCubic.hpp"
#include "VerletListInteractionTemplate.hpp"

namespace espressopp {
namespace interaction {

LOG4ESPP_LOGGER(MultiTabulated::theLogger, "MultiTabulated");

void MultiTabulated::registerTableOnRange(
    const char *filename,
    int itype,
    shared_ptr<analysis::ChemicalConversion> chm,
    real min_value,
    real max_value,
    bool default_table) {


  LOG4ESPP_INFO(theLogger, "registerTable filename=" << filename << " type=" << itype
                            << " min_val=" << min_value << " max_val=" << max_value);
  boost::mpi::communicator world;
  shared_ptr<Interpolation> table;

  switch (itype) {
    case 1:
      table = make_shared<InterpolationLinear>();
      table->read(world, filename);
      break;
    case 2:
      table = make_shared<InterpolationAkima>();
      table->read(world, filename);
      break;
    case 3:
      table = make_shared<InterpolationCubic>();
      table->read(world, filename);
      break;
    default:
      throw std::runtime_error("wrong itype");
  }

  if (default_table) {
    current_table = table;

  }
  chm->onValue.connect(
      boost::bind(&MultiTabulated::onValue, this, min_value, max_value, table, _1)
  );
}

void MultiTabulated::onValue(real min_value, real max_value, shared_ptr<Interpolation> table, real value) {
  if (value >= min_value && value < max_value) {
    current_table = table;

  }
}

/** Register Python interface. */
typedef class VerletListInteractionTemplate<MultiTabulated> VerletListMultiTabulated;

void MultiTabulated::registerPython() {
  using namespace espressopp::python;

  class_<MultiTabulated, bases<Potential> >
      ("interaction_MultiTabulated", init<real>())
      .def("register_table", &MultiTabulated::registerTableOnRange);

  class_<VerletListMultiTabulated, bases<Interaction> >
      ("interaction_VerletListMultiTabulated", init<shared_ptr<VerletList> >())
      .def("setPotential", &VerletListMultiTabulated::setPotential)
      .def("getPotential", &VerletListMultiTabulated::getPotentialPtr);

}

}
}
