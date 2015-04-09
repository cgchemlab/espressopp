/*
  Copyright (C) 2015
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

#include <boost/bind.hpp>
#include <algorithm>
#include <utility>
#include <vector>
#include "storage/Storage.hpp"
#include "FixedPairListType.hpp"
#include "Buffer.hpp"
#include "esutil/Error.hpp"

namespace espressopp {

LOG4ESPP_LOGGER(FixedPairListType::theLogger, "FixedPairListType");

FixedPairListType::FixedPairListType(shared_ptr< storage::Storage > _storage,
    longint type1, longint type2) : FixedPairList(_storage), type1_(type1), type2_(type2) {
  LOG4ESPP_INFO(theLogger, "construct FixedPairListType");
}

bool FixedPairListType::add(longint pid1, longint pid2) {
  if (pid1 > pid2)
    std::swap(pid1, pid2);

  bool returnVal = true;
  System& system = storage->getSystemRef();
  esutil::Error err(system.comm);

  // Particle
  Particle *p1 = storage->lookupRealParticle(pid1);
  Particle *p2 = storage->lookupLocalParticle(pid2);
  if (!p1) {
    // Particle does not exist here, return false
    returnVal = false;
  } else {
    if (!p2) {
      std::stringstream msg;
      msg << "Atomistic bond particle p2 (id=" << pid2 << ") does not exists here "
          << "and cannot be added. "
          << " The pair " << pid1 << " - " << pid2 << " could not be created.";
      err.setException(msg.str());
    }
  }
  err.checkException();

  // Check the type of particles.
  returnVal = returnVal && check_type(p1->type(), p2->type());

  if (returnVal) {
    this->add(p1, p2);
    // ADD THE GLOBAL PAIR
    // see whether the particle already has pairs
    std::pair<GlobalPairs::const_iterator,
      GlobalPairs::const_iterator> equalRange
      = globalPairs.equal_range(pid1);
    if (equalRange.first == globalPairs.end()) {
      // if it hasn't, insert the new pair
      globalPairs.insert(std::make_pair(pid1, pid2));
    } else {
      globalPairs.insert(equalRange.first, std::make_pair(pid1, pid2));
    }
    LOG4ESPP_INFO(theLogger, "added fixed pair to global pair list");
  }
  return returnVal;
}


void FixedPairListType::onParticlesChanged() {
  LOG4ESPP_INFO(theLogger, "rebuild local bond list from global");

  System& system = storage->getSystemRef();
  esutil::Error err(system.comm);

  this->clear();
  longint lastpid1 = -1;
  Particle *p1;
  Particle *p2;

  for (GlobalPairs::const_iterator it = globalPairs.begin(); it != globalPairs.end(); ++it) {
    if (it->first != lastpid1) {
        p1 = storage->lookupRealParticle(it->first);
        if (p1 == NULL) {
          std::stringstream msg;
          msg << "FixedPairListType ";
          msg << "Atomistic bond particle p1 (id=" << it->first << ") does not exists here.";
          err.setException(msg.str());
        }
        lastpid1 = it->first;
    }

    p2 = storage->lookupLocalParticle(it->second);
    if (p2 == NULL) {
      std::stringstream msg;
      msg << "FixedPairListType ";
      msg << "Atomistic bond particle p2 (id=" << it->second << ") does not exists here.";
      err.setException(msg.str());
    }

    // Check types, if valid then add.
    if (check_type(p1->type(), p2->type())) {
      this->add(p1, p2);
    }
  }
  err.checkException();
  LOG4ESPP_INFO(theLogger, "Regenerated local fixed pair list from global list");
}

void FixedPairListType::registerPython() {
  using namespace espressopp::python;  //NOLINT

  bool (FixedPairListType::*pyAdd)(longint pid1, longint pid2)
    = &FixedPairListType::add;

  class_<FixedPairListType, bases<FixedPairList>, shared_ptr<FixedPairListType> >
    ("FixedPairListType", init <shared_ptr<storage::Storage>, int, int >())
      .def("add", pyAdd)
      .def("getBonds",  &FixedPairList::getBonds);
}

}  // end namespace espressopp
