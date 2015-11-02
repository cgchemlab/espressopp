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
#include "NPart.hpp"
#include "storage/DomainDecomposition.hpp"
#include <boost/unordered_set.hpp>
#include "iterator/CellListIterator.hpp"

using namespace espressopp;

namespace espressopp {
  namespace analysis {
    using namespace iterator;
    real NPart::compute_real() const {

      int myN, systemN;
      System& system = getSystemRef();
      if (!has_types) {
        myN = system.storage->getNRealParticles();
        myN += system.storage->getNAdressParticles();
      } else {
        myN = 0;
        CellList realCells = system.storage->getRealCells();
        bool has_adress = (bool) system.storage->getFixedTuples();
        if (has_adress) {
          FixedTupleListAdress::iterator it2;
          shared_ptr<FixedTupleListAdress> fixedtupleList = system.storage->getFixedTuples();
          std::vector<Particle*> atList;
          for (CellListIterator cit(realCells); !cit.isDone(); ++cit) {
            if (valid_type_ids.count(cit->type()))
              myN++;
            Particle &vp = *cit;
            it2 = fixedtupleList->find(&vp);
            if (it2 != fixedtupleList->end()) {
              atList = it2->second;
              for (std::vector<Particle*>::iterator it3 = atList.begin(); it3 != atList.end(); ++it3) {
                Particle &at = **it3;
                if (valid_type_ids.count(at.type()))
                  myN++;
              }
            }
          }
        } else {
          for (CellListIterator cit(realCells); !cit.isDone(); ++cit) {
            if (valid_type_ids.count(cit->type()))
              myN++;
          }
        }
      }
      boost::mpi::reduce(*getSystem()->comm, myN, systemN, std::plus<int>(), 0);
      
      return 1.0*systemN;
    }

    void NPart::addType(longint type_id) {
      has_types = true;
      valid_type_ids.insert(type_id);
    }

    bool NPart::removeType(longint type_id) {
      bool ret_val = valid_type_ids.erase(type_id);
      has_types = valid_type_ids.size() > 0;
      return ret_val;
    }

    void NPart::registerPython() {
      using namespace espressopp::python;
      class_<NPart, bases< Observable > >
        ("analysis_NPart", init< shared_ptr< System > >())
        .def("add_type", &NPart::addType)
        .def("remove_type", &NPart::removeType)
      ;
    }
  }
}
