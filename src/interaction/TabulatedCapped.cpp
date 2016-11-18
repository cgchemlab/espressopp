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
#include "TabulatedCapped.hpp"
#include "InterpolationLinear.hpp"
#include "InterpolationAkima.hpp"
#include "InterpolationCubic.hpp"
#include "VerletListInteractionTemplate.hpp"
#include "VerletListAdressInteractionTemplate.hpp"
#include "VerletListHadressInteractionTemplate.hpp"
#include "VerletListDynamicResolutionInteractionTemplate.hpp"
#include "CellListAllPairsInteractionTemplate.hpp"
#include "FixedPairListInteractionTemplate.hpp"
#include "FixedPairListTypesInteractionTemplate.hpp"
#include "FixedPairListLambdaInteractionTemplate.hpp"
#include "FixedPairListTypesLambdaInteractionTemplate.hpp"

namespace espressopp {
namespace interaction {


void TabulatedCapped::setFilename(int itype, const char *_filename) {
  boost::mpi::communicator world;
  filename = _filename;

  if (itype == 1) { // create a new InterpolationLinear
    table = make_shared<InterpolationLinear>();
    table->read(world, _filename);
  } else if (itype == 2) { // create a new InterpolationAkima
    table = make_shared<InterpolationAkima>();
    table->read(world, _filename);
  } else if (itype == 3) { // create a new InterpolationCubic
    table = make_shared<InterpolationCubic>();
    table->read(world, _filename);
  }
}

real TabulatedCapped::getCaprad() const { return caprad_; }

typedef class VerletListInteractionTemplate<TabulatedCapped> VerletListTabulatedCapped;
typedef class VerletListAdressInteractionTemplate<TabulatedCapped, TabulatedCapped> VerletListAdressTabulatedCapped;
typedef class VerletListHadressInteractionTemplate<TabulatedCapped, TabulatedCapped> VerletListHadressTabulatedCapped;
typedef class VerletListDynamicResolutionInteractionTemplate<TabulatedCapped>
    VerletListDynamicResolutionTabulatedCapped;
typedef class CellListAllPairsInteractionTemplate<TabulatedCapped> CellListTabulatedCapped;
typedef class FixedPairListInteractionTemplate<TabulatedCapped> FixedPairListTabulatedCapped;
typedef class FixedPairListTypesInteractionTemplate<TabulatedCapped> FixedPairListTypesTabulatedCapped;
typedef class FixedPairListLambdaInteractionTemplate<TabulatedCapped> FixedPairListLambdaTabulatedCapped;
typedef class FixedPairListTypesLambdaInteractionTemplate<TabulatedCapped> FixedPairListTypesLambdaTabulatedCapped;

//////////////////////////////////////////////////
// REGISTRATION WITH PYTHON
//////////////////////////////////////////////////
void TabulatedCapped::registerPython() {
  using namespace espressopp::python;

  class_<TabulatedCapped, bases<Potential> >
      ("interaction_TabulatedCapped", init<int, const char *, real, real>())
      .def(init<int, const char *>())
      .def(init<int, const char *, real>())
      .add_property("filename", &TabulatedCapped::getFilename, &TabulatedCapped::setFilename)
      .add_property("itype", make_getter(&TabulatedCapped::interpolationType))
      .add_property("caprad", make_getter(&TabulatedCapped::caprad_))
      .def_pickle(TabulatedCapped_pickle());

  class_<VerletListTabulatedCapped, bases<Interaction> >
      ("interaction_VerletListTabulatedCapped", init<shared_ptr<VerletList> >())
      .def("setPotential", &VerletListTabulatedCapped::setPotential)
      .def("getPotential", &VerletListTabulatedCapped::getPotentialPtr);

  class_<VerletListAdressTabulatedCapped, bases<Interaction> >
      ("interaction_VerletListAdressTabulatedCapped",
       init<shared_ptr<VerletListAdress>,
            shared_ptr<FixedTupleListAdress> >()
      )
      .def("setPotentialAT", &VerletListAdressTabulatedCapped::setPotentialAT)
      .def("setPotentialCG", &VerletListAdressTabulatedCapped::setPotentialCG);;

  class_<VerletListHadressTabulatedCapped, bases<Interaction> >
      ("interaction_VerletListHadressTabulatedCapped",
       init<shared_ptr<VerletListAdress>,
            shared_ptr<FixedTupleListAdress> >()
      )
      .def("setPotentialAT", &VerletListHadressTabulatedCapped::setPotentialAT)
      .def("setPotentialCG", &VerletListHadressTabulatedCapped::setPotentialCG);;

  class_<VerletListDynamicResolutionTabulatedCapped, bases<Interaction> >
      ("interaction_VerletListDynamicResolutionTabulatedCapped", init<shared_ptr<VerletList>, bool>())
      .def("getVerletList", &VerletListDynamicResolutionTabulatedCapped::getVerletList)
      .def("setPotential", &VerletListDynamicResolutionTabulatedCapped::setPotential)
      .def("getPotential", &VerletListDynamicResolutionTabulatedCapped::getPotentialPtr)
      .def("setMaxForce", &VerletListDynamicResolutionTabulatedCapped::setMaxForce);

  class_<CellListTabulatedCapped, bases<Interaction> >
      ("interaction_CellListTabulatedCapped", init<shared_ptr<storage::Storage> >())
      .def("setPotential", &CellListTabulatedCapped::setPotential);;

  class_<FixedPairListTabulatedCapped, bases<Interaction> >
      ("interaction_FixedPairListTabulatedCapped",
       init<shared_ptr<System>,
            shared_ptr<FixedPairList>,
            shared_ptr<TabulatedCapped> >()
      )
      .def(init<shared_ptr<System>, shared_ptr<FixedPairListAdress>, shared_ptr<TabulatedCapped> >())
      .def("setPotential", &FixedPairListTabulatedCapped::setPotential)
      .def("getPotential", &FixedPairListTabulatedCapped::getPotential)
      .def("setFixedPairList", &FixedPairListTabulatedCapped::setFixedPairList)
      .def("getFixedPairList", &FixedPairListTabulatedCapped::getFixedPairList);;

  class_<FixedPairListTypesTabulatedCapped, bases<Interaction> >
      ("interaction_FixedPairListTypesTabulatedCapped",
       init<shared_ptr<System>, shared_ptr<FixedPairList> >())
      .def(init<shared_ptr<System>, shared_ptr<FixedPairListAdress> >())
      .def("setPotential", &FixedPairListTypesTabulatedCapped::setPotential)
      .def("getPotential", &FixedPairListTypesTabulatedCapped::getPotentialPtr)
      .def("setFixedPairList", &FixedPairListTypesTabulatedCapped::setFixedPairList)
      .def("getFixedPairList", &FixedPairListTypesTabulatedCapped::getFixedPairList);

  class_<FixedPairListLambdaTabulatedCapped, bases<Interaction> >
      ("interaction_FixedPairListLambdaTabulatedCapped",
       init<shared_ptr<System>,
            shared_ptr<FixedPairListLambda>,
            shared_ptr<TabulatedCapped> >()
      )
      .def("setPotential", &FixedPairListLambdaTabulatedCapped::setPotential)
      .def("getPotential", &FixedPairListLambdaTabulatedCapped::getPotential)
      .def("setFixedPairList", &FixedPairListLambdaTabulatedCapped::setFixedPairList)
      .def("getFixedPairList", &FixedPairListLambdaTabulatedCapped::getFixedPairList);;

  class_<FixedPairListTypesLambdaTabulatedCapped, bases<Interaction> >
      ("interaction_FixedPairListTypesLambdaTabulatedCapped",
       init<shared_ptr<System>, shared_ptr<FixedPairListLambda> >())
      .def("setPotential", &FixedPairListTypesLambdaTabulatedCapped::setPotential)
      .def("getPotential", &FixedPairListTypesLambdaTabulatedCapped::getPotentialPtr)
      .def("setFixedPairList", &FixedPairListTypesLambdaTabulatedCapped::setFixedPairList)
      .def("getFixedPairList", &FixedPairListTypesLambdaTabulatedCapped::getFixedPairList);
}

}
}
