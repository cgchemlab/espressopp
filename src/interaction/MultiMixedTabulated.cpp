/*
  Copyright (C) 2017
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
#include "MultiMixedTabulated.hpp"
#include "InterpolationLinear.hpp"
#include "InterpolationAkima.hpp"
#include "InterpolationCubic.hpp"
#include "VerletListInteractionTemplate.hpp"

namespace espressopp {
namespace interaction {

LOG4ESPP_LOGGER(MultiMixedTabulated::theLogger, "MultiMixedTabulated");

void MultiMixedTabulated::registerTableOnRange(
    std::string tab1,
    std::string tab2,
    shared_ptr<analysis::ChemicalConversion> chm,
    real min_value,
    real max_value) {

  LOG4ESPP_INFO(theLogger,
                "registerTable tab1=" << tab1 << " tab2=" << tab2 << " type=" << interpolationType
                << " min_val=" << min_value << " max_val=" << max_value);
  boost::mpi::communicator world;
  shared_ptr<Interpolation> table1;
  shared_ptr<Interpolation> table2;

  switch (interpolationType) {
    case 1:
      table1 = make_shared<InterpolationLinear>();
      table1->read(world, tab1.c_str());
      table2 = make_shared<InterpolationLinear>();
      table2->read(world, tab2.c_str());
      break;
    case 2:
      table1 = make_shared<InterpolationAkima>();
      table1->read(world, tab1.c_str());
      table2 = make_shared<InterpolationAkima>();
      table2->read(world, tab2.c_str());
      break;
    case 3:
      table1 = make_shared<InterpolationCubic>();
      table1->read(world, tab1.c_str());
      table2 = make_shared<InterpolationCubic>();
      table2->read(world, tab2.c_str());
      break;
    default:
      throw std::runtime_error("wrong itype");
  }

  if (!current_table1)
    current_table1 = table1;

  if (!current_table2)
    current_table2 = table2;

  chm->onValue.connect(
      boost::bind(&MultiMixedTabulated::onValue, this, min_value, max_value, table1, table2, _1)
  );
}

void MultiMixedTabulated::onValue(real min_value, real max_value, shared_ptr<Interpolation> table1, shared_ptr<Interpolation> table2, real value) {
  if (value > min_value && value <= max_value) {
    current_table1 = table1;
    current_table2 = table2;
  }
  mix_value_ = value;
}

/** Register Python interface. */
typedef class VerletListInteractionTemplate<MultiMixedTabulated> VerletListMultiMixedTabulated;

void MultiMixedTabulated::registerPython() {
  using namespace espressopp::python;

  class_<MultiMixedTabulated, bases<Potential> >
      ("interaction_MultiMixedTabulated", init<int, real>())
      .def("register_table", &MultiMixedTabulated::registerTableOnRange);

  class_<VerletListMultiMixedTabulated, bases<Interaction> >
      ("interaction_VerletListMultiMixedTabulated", init<shared_ptr<VerletList> >())
      .def("setPotential", &VerletListMultiMixedTabulated::setPotential)
      .def("getPotential", &VerletListMultiMixedTabulated::getPotentialPtr);

}

}
}
