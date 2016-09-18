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
#include "MixedTabulated.hpp"
#include "VerletListInteractionTemplate.hpp"

namespace espressopp {
namespace interaction {

LOG4ESPP_LOGGER(MixedTabulated::theLogger, "MixedTabulated");
LOG4ESPP_LOGGER(MixedTabulated::theLocalLogger, "MixedTabulatedLocal");


void MixedTabulated::onValue(real value) {
  LOG4ESPP_DEBUG(theLocalLogger, "previous mixed value=" << mix_value_ << " new=" << value);
  mix_value_ = value;
}

/** Register Python interface. */
typedef class VerletListInteractionTemplate<MixedTabulated> VerletListMixedTabulated;

void MixedTabulated::registerPython() {
  using namespace espressopp::python;

  class_<MixedTabulated, bases<Potential> >
      ("interaction_MixedTabulated",
          init<int, const char*, const char*, shared_ptr<analysis::ChemicalConversion>, real, real>())
          .def(init<int, const char*, const char*, real, real>())  // static mixing
          .add_property("mix_value",
                        &MixedTabulated::mix_value,
                        &MixedTabulated::set_mix_value);

  class_<VerletListMixedTabulated, bases<Interaction> >
      ("interaction_VerletListMixedTabulated", init<shared_ptr<VerletList> >())
      .def("setPotential", &VerletListMixedTabulated::setPotential)
      .def("getPotential", &VerletListMixedTabulated::getPotentialPtr);
}

}
}
