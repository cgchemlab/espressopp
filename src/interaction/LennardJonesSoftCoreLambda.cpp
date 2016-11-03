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
#include "LennardJonesSoftCoreLambda.hpp"
#include "Tabulated.hpp"
#include "VerletListInteractionTemplate.hpp"

namespace espressopp {
namespace interaction {

typedef class VerletListInteractionTemplate<LennardJonesSoftCoreLambda> VerletListLennardJonesSoftCoreLambda;

LOG4ESPP_LOGGER(LennardJonesSoftCoreLambda::theLogger, "LennardJonesSoftCoreLambda");

void LennardJonesSoftCoreLambda::registerPython() {
  using namespace espressopp::python;

  class_<LennardJonesSoftCoreLambda, bases<Potential> >
      ("interaction_LennardJonesSoftCoreLambda", init<real, real, real, real, real, real>())
      .def(init<real, real, real, real, real>())
      .add_property("sigma", &LennardJonesSoftCoreLambda::getSigma, &LennardJonesSoftCoreLambda::setSigma)
      .add_property("epsilon", &LennardJonesSoftCoreLambda::getEpsilon, &LennardJonesSoftCoreLambda::setEpsilon)
      .add_property("alpha", &LennardJonesSoftCoreLambda::getAlpha, &LennardJonesSoftCoreLambda::setAlpha)
      .add_property("weight_power", &LennardJonesSoftCoreLambda::getWeightPower,
                    &LennardJonesSoftCoreLambda::setWeightPower)
      .def_pickle(LennardJonesSoftCoreLambda_pickle());

  class_<VerletListLennardJonesSoftCoreLambda, bases<Interaction> >
      ("interaction_VerletListLennardJonesSoftCoreLambda", init<shared_ptr<VerletList> >())
      .def("getVerletList", &VerletListLennardJonesSoftCoreLambda::getVerletList)
      .def("setPotential", &VerletListLennardJonesSoftCoreLambda::setPotential)
      .def("getPotential", &VerletListLennardJonesSoftCoreLambda::getPotentialPtr);
}

}
}
