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
#include "LennardJonesForceCapped.hpp"
#include "Tabulated.hpp"
#include "VerletListInteractionTemplate.hpp"
#include "VerletListAdressInteractionTemplate.hpp"
#include "VerletListHadressInteractionTemplate.hpp"
#include "CellListAllPairsInteractionTemplate.hpp"
#include "FixedPairListInteractionTemplate.hpp"

namespace espressopp {
  namespace interaction {

    typedef class VerletListInteractionTemplate <LennardJonesForceCapped>
        VerletListLennardJonesForceCapped;
    typedef class VerletListAdressInteractionTemplate <LennardJonesForceCapped, Tabulated>
        VerletListAdressLennardJonesForceCapped;
    typedef class VerletListHadressInteractionTemplate <LennardJonesForceCapped, Tabulated>
        VerletListHadressLennardJonesForceCapped;
    typedef class CellListAllPairsInteractionTemplate <LennardJonesForceCapped>
        CellListLennardJonesForceCapped;
    typedef class FixedPairListInteractionTemplate <LennardJonesForceCapped>
        FixedPairListLennardJonesForceCapped;

    //////////////////////////////////////////////////
    // REGISTRATION WITH PYTHON
    //////////////////////////////////////////////////
    void 
    LennardJonesForceCapped::registerPython() {
      using namespace espressopp::python;

      class_< LennardJonesForceCapped, bases< Potential > >
    	("interaction_LennardJonesForceCapped", init< real, real, real, real >())
	    .def(init< real, real, real, real, real>())
    	.add_property("sigma", &LennardJonesForceCapped::getSigma, &LennardJonesForceCapped::setSigma)
    	.add_property("epsilon", &LennardJonesForceCapped::getEpsilon, &LennardJonesForceCapped::setEpsilon)
    	.add_property("max_force", &LennardJonesForceCapped::getMaxforce, &LennardJonesForceCapped::setMaxforce)
    	.def_pickle(LennardJonesForceCapped_pickle())
      ;

      class_< VerletListLennardJonesForceCapped, bases< Interaction > >
        ("interaction_VerletListLennardJonesForceCapped", init< shared_ptr<VerletList> >())
        .def("setPotential", &VerletListLennardJonesForceCapped::setPotential, return_value_policy< reference_existing_object >())
        .def("getPotential", &VerletListLennardJonesForceCapped::getPotential, return_value_policy< reference_existing_object >())
      ;

      class_< VerletListAdressLennardJonesForceCapped, bases< Interaction > >
        ("interaction_VerletListAdressLennardJonesForceCapped",
         init< shared_ptr<VerletListAdress>, shared_ptr<FixedTupleListAdress> >())
         .def("setPotentialAT", &VerletListAdressLennardJonesForceCapped::setPotentialAT)
         .def("setPotentialCG", &VerletListAdressLennardJonesForceCapped::setPotentialCG)
         .def("getPotentialAT", &VerletListAdressLennardJonesForceCapped::getPotentialAT,
                 return_value_policy< reference_existing_object >())
         .def("getPotentialCG", &VerletListAdressLennardJonesForceCapped::getPotentialCG,
                       return_value_policy< reference_existing_object >());
      ;

      class_< VerletListHadressLennardJonesForceCapped, bases< Interaction > >
        ("interaction_VerletListHadressLennardJonesForceCapped",
         init< shared_ptr<VerletListAdress>, shared_ptr<FixedTupleListAdress> >())
         .def("setPotentialAT", &VerletListHadressLennardJonesForceCapped::setPotentialAT)
         .def("setPotentialCG", &VerletListHadressLennardJonesForceCapped::setPotentialCG)
         .def("getPotentialAT", &VerletListHadressLennardJonesForceCapped::getPotentialAT,
                 return_value_policy< reference_existing_object >())
         .def("getPotentialCG", &VerletListHadressLennardJonesForceCapped::getPotentialCG,
                       return_value_policy< reference_existing_object >());
      ;
      
      class_< CellListLennardJonesForceCapped, bases< Interaction > >
        ("interaction_CellListLennardJonesForceCapped", init< shared_ptr< storage::Storage > >())
        .def("setPotential", &CellListLennardJonesForceCapped::setPotential)
        .def("getPotential", &CellListLennardJonesForceCapped::getPotential, return_value_policy< reference_existing_object >());
	  ;

      class_< FixedPairListLennardJonesForceCapped, bases< Interaction > >
        ("interaction_FixedPairListLennardJonesForceCapped",
          init< shared_ptr<System>, shared_ptr<FixedPairList>, shared_ptr<LennardJonesForceCapped> >())
          .def(init< shared_ptr<System>, shared_ptr<FixedPairListAdress>, shared_ptr<LennardJonesForceCapped> >())
          .def("setPotential", &FixedPairListLennardJonesForceCapped::setPotential);
      ;
    }
    
  }
}
