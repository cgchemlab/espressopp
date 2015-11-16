/*
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
#include "LennardJonesLambda.hpp"
#include "Tabulated.hpp"
#include "VerletListInteractionTemplate.hpp"
#include "VerletListNonReciprocalInteractionTemplate.hpp"
#include "VerletListAdressInteractionTemplate.hpp"
#include "VerletListHadressInteractionTemplate.hpp"
#include "VerletListDynamicResolutionInteractionTemplate.hpp"
#include "CellListAllPairsInteractionTemplate.hpp"
#include "FixedPairListInteractionTemplate.hpp"

namespace espressopp {
  namespace interaction {

    typedef class VerletListInteractionTemplate <LennardJonesLambda>
        VerletListLennardJonesLambda;

    LOG4ESPP_LOGGER(LennardJonesLambda::theLogger, "LennardJonesLambda");

    //////////////////////////////////////////////////
    // REGISTRATION WITH PYTHON
    //////////////////////////////////////////////////
    void
    LennardJonesLambda::registerPython() {
      using namespace espressopp::python;

      class_< LennardJonesLambda, bases< Potential > >
    	("interaction_LennardJonesLambda", init< real, real, real >())
	    .def(init< real, real, real, real >())
    	.add_property("sigma", &LennardJonesLambda::getSigma, &LennardJonesLambda::setSigma)
    	.add_property("epsilon", &LennardJonesLambda::getEpsilon, &LennardJonesLambda::setEpsilon)
        .def_pickle(LennardJonesLambda_pickle())
      ;

      class_< VerletListLennardJonesLambda, bases< Interaction > >
        ("interaction_VerletListLennardJonesLambda", init< shared_ptr<VerletList> >())
        .def("getVerletList", &VerletListLennardJonesLambda::getVerletList)
        .def("setPotential", &VerletListLennardJonesLambda::setPotential)
        .def("getPotential", &VerletListLennardJonesLambda::getPotentialPtr)
      ;
    }

  }
}
