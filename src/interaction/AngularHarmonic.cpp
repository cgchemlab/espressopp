/*
  Copyright (C) 2016
      Jakub Krajniak (jkrajniak at gmail.com)
  Copyright (C) 2012,2013
      Max Planck Institute for Polymer Research
  Copyright (C) 2008,2009,2010,2011
      Max-Planck-Institute for Polymer Research & Fraunhofer SCAI
  
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
#include "AngularHarmonic.hpp"
#include <FixedTripleListLambda.hpp>
#include "FixedTripleListInteractionTemplate.hpp"
#include "FixedTripleListTypesInteractionTemplate.hpp"
#include "FixedTripleListLambdaInteractionTemplate.hpp"
#include "FixedTripleListTypesLambdaInteractionTemplate.hpp"

namespace espressopp {
  namespace interaction {
    //////////////////////////////////////////////////
    // REGISTRATION WITH PYTHON
    //////////////////////////////////////////////////
    void 
    AngularHarmonic::registerPython() {
      using namespace espressopp::python;

      class_< AngularHarmonic, bases< AngularPotential > >
    	("interaction_AngularHarmonic", init< real, real >())
	.add_property("K", &AngularHarmonic::getK, &AngularHarmonic::setK)
	.add_property("theta0", &AngularHarmonic::getTheta0, &AngularHarmonic::setTheta0)
    	;

    typedef class FixedTripleListInteractionTemplate<AngularHarmonic>
        FixedTripleListAngularHarmonic;

    typedef class FixedTripleListTypesInteractionTemplate<AngularHarmonic>
        FixedTripleListTypesAngularHarmonic;

    typedef class FixedTripleListLambdaInteractionTemplate<AngularHarmonic>
        FixedTripleListLambdaAngularHarmonic;

    typedef class FixedTripleListTypesLambdaInteractionTemplate<AngularHarmonic>
        FixedTripleListTypesLambdaAngularHarmonic;
        
      class_ <FixedTripleListAngularHarmonic, bases <Interaction> >
        ("interaction_FixedTripleListAngularHarmonic",
           init<shared_ptr<System>,
                shared_ptr<FixedTripleList>,
                shared_ptr<AngularHarmonic> >())
        .def(init< shared_ptr<System>, shared_ptr<FixedTripleListAdress>, shared_ptr<AngularHarmonic> >())
        .def("setPotential", &FixedTripleListAngularHarmonic::setPotential)
        .def("getFixedTripleList", &FixedTripleListAngularHarmonic::getFixedTripleList);

      class_ <FixedTripleListLambdaAngularHarmonic, bases <Interaction> >
          ("interaction_FixedTripleListLambdaAngularHarmonic",
           init<shared_ptr<System>,
                shared_ptr<FixedTripleListLambda>,
                shared_ptr<AngularHarmonic> >())
          .def("setPotential", &FixedTripleListLambdaAngularHarmonic::setPotential)
          .def("getFixedTripleList", &FixedTripleListLambdaAngularHarmonic::getFixedTripleList);

      class_< FixedTripleListTypesAngularHarmonic, bases< Interaction > >
          ("interaction_FixedTripleListTypesAngularHarmonic",
           init< shared_ptr<System>, shared_ptr<FixedTripleList> >())
          .def("setPotential", &FixedTripleListTypesAngularHarmonic::setPotential)
          .def("getPotential", &FixedTripleListTypesAngularHarmonic::getPotentialPtr)
          .def("setFixedTripleList", &FixedTripleListTypesAngularHarmonic::setFixedTripleList)
          .def("getFixedTripleList", &FixedTripleListTypesAngularHarmonic::getFixedTripleList);

      class_< FixedTripleListTypesLambdaAngularHarmonic, bases< Interaction > >
          ("interaction_FixedTripleListTypesLambdaAngularHarmonic",
           init< shared_ptr<System>, shared_ptr<FixedTripleListLambda> >())
          .def("setPotential", &FixedTripleListTypesLambdaAngularHarmonic::setPotential)
          .def("getPotential", &FixedTripleListTypesLambdaAngularHarmonic::getPotentialPtr)
          .def("setFixedTripleList", &FixedTripleListTypesLambdaAngularHarmonic::setFixedTripleList)
          .def("getFixedTripleList", &FixedTripleListTypesLambdaAngularHarmonic::getFixedTripleList);
    }
  }
}
