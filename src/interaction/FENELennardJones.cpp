/*
  Copyright (C) 2017
      Jakub Krajniak (jkrajnika at gmail.com)
  
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
#include "FENELennardJones.hpp"
#include "FixedPairListInteractionTemplate.hpp"
#include "FixedPairListLambdaInteractionTemplate.hpp"
#include "FixedPairListTypesLambdaInteractionTemplate.hpp"
#include "FixedPairListTypesInteractionTemplate.hpp"


namespace espressopp {
  namespace interaction {
    typedef class FixedPairListInteractionTemplate< FENELennardJones > FixedPairListFENELennardJones;

    typedef class FixedPairListLambdaInteractionTemplate<FENELennardJones>
      FixedPairListLambdaFENELennardJones;
    typedef class FixedPairListLambdaInteractionTemplate<FENELennardJones>
        FixedPairListLambdaFENELennardJones;
    typedef class FixedPairListTypesLambdaInteractionTemplate<FENELennardJones>
        FixedPairListTypesLambdaFENELennardJones;
    typedef class FixedPairListTypesInteractionTemplate<FENELennardJones>
      FixedPairListTypesFENELennardJones;

    LOG4ESPP_LOGGER(FENELennardJones::theLogger, "FENELennardJones");
    //////////////////////////////////////////////////
    // REGISTRATION WITH PYTHON
    //////////////////////////////////////////////////
    void 
    FENELennardJones::registerPython() {
      using namespace espressopp::python;
      class_< FENELennardJones, bases< Potential > >
      ("interaction_FENELennardJones", init< real, real, real, real, real, real >())
	  .def(init< real, real, real, real, real, real, real >())
	  .add_property("K", &FENELennardJones::getK, &FENELennardJones::setK)
	  .add_property("r0", &FENELennardJones::getR0, &FENELennardJones::setR0)
	  .add_property("rMax", &FENELennardJones::getRMax, &FENELennardJones::setRMax)
      .add_property("sigma", &FENELennardJones::getSigma, &FENELennardJones::setSigma)
      .add_property("epsilon", &FENELennardJones::getEpsilon, &FENELennardJones::setEpsilon)
	  ;

      class_< FixedPairListFENELennardJones, bases< Interaction > >
      ("interaction_FixedPairListFENELennardJones", init< shared_ptr<System>, shared_ptr<FixedPairList>, shared_ptr<FENELennardJones> >())
      .def(init< shared_ptr<System>, shared_ptr<FixedPairListAdress>, shared_ptr<FENELennardJones> >())
      .def("setPotential", &FixedPairListFENELennardJones::setPotential)
      .def("getPotential", &FixedPairListFENELennardJones::getPotential)
      .def("setFixedPairList", &FixedPairListFENELennardJones::setFixedPairList)
      .def("getFixedPairList", &FixedPairListFENELennardJones::getFixedPairList)
      ;

      class_< FixedPairListTypesFENELennardJones, bases< Interaction > >
          ("interaction_FixedPairListTypesFENELennardJones",
           init< shared_ptr<System>, shared_ptr<FixedPairList> >())
          .def(init< shared_ptr<System>, shared_ptr<FixedPairListAdress> >())
          .def("setPotential", &FixedPairListTypesFENELennardJones::setPotential)
          .def("getPotential", &FixedPairListTypesFENELennardJones::getPotentialPtr)
          .def("setFixedPairList", &FixedPairListTypesFENELennardJones::setFixedPairList)
          .def("getFixedPairList", &FixedPairListTypesFENELennardJones::getFixedPairList);
      ;
      class_< FixedPairListLambdaFENELennardJones, bases< Interaction > >
          ("interaction_FixedPairListLambdaFENELennardJones",
           init< shared_ptr<System>, shared_ptr<FixedPairListLambda>, shared_ptr<FENELennardJones> >())
          .def("setPotential", &FixedPairListLambdaFENELennardJones::setPotential)
          .def("getPotential", &FixedPairListLambdaFENELennardJones::getPotential)
          .def("setFixedPairList", &FixedPairListLambdaFENELennardJones::setFixedPairList)
          .def("setFixedPairList", &FixedPairListLambdaFENELennardJones::setFixedPairList);
      ;

      class_< FixedPairListTypesLambdaFENELennardJones, bases< Interaction > >
          ("interaction_FixedPairListTypesLambdaFENELennardJones",
           init< shared_ptr<System>, shared_ptr<FixedPairListLambda> >())
          .def("setPotential", &FixedPairListTypesLambdaFENELennardJones::setPotential)
          .def("getPotential", &FixedPairListTypesLambdaFENELennardJones::getPotentialPtr)
          .def("setFixedPairList", &FixedPairListTypesLambdaFENELennardJones::setFixedPairList)
          .def("getFixedPairList", &FixedPairListTypesLambdaFENELennardJones::getFixedPairList);
      ;

    }

  }
}
