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
#include "LennardJonesUnique.hpp"
#include "Tabulated.hpp"
#include "VerletListInteractionTemplate.hpp"

namespace espressopp {
  namespace interaction {

    typedef class VerletListInteractionTemplate <LennardJonesUnique>
        VerletListLennardJonesUnique;
    LOG4ESPP_LOGGER(LennardJonesUnique::theLogger, "LennardJonesUnique");

    //////////////////////////////////////////////////
    // REGISTRATION WITH PYTHON
    //////////////////////////////////////////////////
    void
    LennardJonesUnique::registerPython() {
      using namespace espressopp::python;

      class_< LennardJonesUnique, bases< Potential > >
    	("interaction_LennardJonesUnique", init< real, real, real >())
	    .def(init< real, real, real, real >())
    	.add_property("sigma", &LennardJonesUnique::getSigma, &LennardJonesUnique::setSigma)
    	.add_property("epsilon", &LennardJonesUnique::getEpsilon, &LennardJonesUnique::setEpsilon)
        .def_pickle(LennardJonesUnique_pickle())
      ;
      
      class_< VerletListLennardJonesUnique, bases< Interaction > >
        ("interaction_VerletListLennardJonesUnique", init< shared_ptr<VerletList> >())
        .def("getVerletList", &VerletListLennardJonesUnique::getVerletList)
        .def("setPotential", &VerletListLennardJonesUnique::setPotential)
        .def("getPotential", &VerletListLennardJonesUnique::getPotentialPtr)
      ;

    }

  }
}
