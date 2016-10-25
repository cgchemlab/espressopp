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
#include "FixedPairLambdaList.hpp"
#include "storage/Storage.hpp"
#include "bc/BC.hpp"
#include "esutil/Error.hpp"

namespace espressopp {

LOG4ESPP_LOGGER(FixedPairLambdaList::theLogger, "FixedPairLambdaList");


FixedPairLambdaList::FixedPairLambdaList(shared_ptr<storage::Storage> _storage, real initLambda)
    : storage(_storage), pairsLambda(), initLambda_(initLambda) {
  LOG4ESPP_INFO(theLogger, "construct FixedPairLambdaList");

  con1 = storage->beforeSendParticles.connect
      (boost::bind(&FixedPairLambdaList::beforeSendParticles, this, _1, _2));
  con2 = storage->afterRecvParticles.connect
      (boost::bind(&FixedPairLambdaList::afterRecvParticles, this, _1, _2));
  con3 = storage->onParticlesChanged.connect
      (boost::bind(&FixedPairLambdaList::onParticlesChanged, this));
}

FixedPairLambdaList::~FixedPairLambdaList() {
  LOG4ESPP_INFO(theLogger, "~FixedPairLambdaList");
  con1.disconnect();
  con2.disconnect();
  con3.disconnect();
}

bool FixedPairLambdaList::add(longint pid1, longint pid2) {
  bool returnVal = true;
  System &system = storage->getSystemRef();
  esutil::Error err(system.comm);

  if (pid1 > pid2)
    std::swap(pid1, pid2);

  // ADD THE LOCAL PAIR
  Particle *p1 = storage->lookupRealParticle(pid1);
  Particle *p2 = storage->lookupLocalParticle(pid2);
  if (!p1) {
    returnVal = false;
  } else {
    if (!p2) {
      std::stringstream msg;
      msg << "bond particle p2 " << pid2 << " does not exists here and cannot be added";
      err.setException(msg.str());
    }
  }
  err.checkException();

  if (returnVal) {
    bool found = false;
    std::pair<PairsLambda::const_iterator,
              PairsLambda::const_iterator> equalRange = pairsLambda.equal_range(pid1);
    if (equalRange.first != pairsLambda.end()) {
      for (PairsLambda::const_iterator it = equalRange.first; it != equalRange.second && !found; ++it) {
        if (it->second.first == pid2)
          found = true;
      }
    }
    returnVal = !found;
    if (!found) {
      this->push_back(ParticlePairLambda(p1, p2, initLambda_));
      onTupleAdded(pid1, pid2);
      pairsLambda.insert(std::make_pair(pid1, std::make_pair(pid2, initLambda_)));
      LOG4ESPP_INFO(theLogger, "added fixed pair to global pair lambda list");
    }
  }
  return returnVal;
}

python::list FixedPairLambdaList::getPairs() {
  python::tuple pair;
  python::list pairs;
  for (PairsLambda::const_iterator it = pairsLambda.begin(); it != pairsLambda.end(); it++) {
    pair = python::make_tuple(it->first, it->second.first);
    pairs.append(pair);
  }

  return pairs;
}

python::list FixedPairLambdaList::getPairsLambda() {
  python::tuple pair;
  python::list pairs;
  for (PairsLambda::const_iterator it = pairsLambda.begin(); it != pairsLambda.end(); it++) {
    pair = python::make_tuple(it->first, it->second.first, it->second.second);
    pairs.append(pair);
  }

  return pairs;
}

real FixedPairLambdaList::getLambda(longint pid1, longint pid2) {
  real returnVal = -3;

  PairsLambda::iterator itr;
  PairsLambda::iterator lastElement;

  // locate an iterator to the first pair object associated with key
  itr = pairsLambda.find(pid1);
  if (itr == pairsLambda.end())
    return returnVal; // no elements associated with key, so return immediately

  // get an iterator to the element that is one past the last element associated with key
  lastElement = pairsLambda.upper_bound(pid1);

  for (; itr != lastElement; ++itr) {
    if (pid2 == itr->second.first) {
      returnVal = itr->second.second;
      break;
    }
  }

  return returnVal;
}

void FixedPairLambdaList::beforeSendParticles(ParticleList &pl, OutBuffer &buf) {
  std::vector<longint> toSendInt;
  std::vector<real> toSendReal;
  // loop over the particle list
  for (ParticleList::Iterator pit(pl); pit.isValid(); ++pit) {
    longint pid = pit->id();

    LOG4ESPP_DEBUG(theLogger, "send particle with pid " << pid << ", find pairs");
    int n = pairsLambda.count(pid);

    if (n > 0) {
      std::pair<PairsLambda::const_iterator,
                PairsLambda::const_iterator> equalRange = pairsLambda.equal_range(pid);

      // first write the pid of the first particle
      // then the number of partners
      // and then the pids of the partners
      toSendInt.reserve(toSendInt.size() + n + 1);
      toSendReal.reserve(toSendReal.size() + n);
      toSendInt.push_back(pid);
      toSendInt.push_back(n);
      for (PairsLambda::const_iterator it = equalRange.first;
           it != equalRange.second; ++it) {
        toSendInt.push_back(it->second.first);
        toSendReal.push_back(it->second.second);
      }
      // delete all of these pairs from the global list
      pairsLambda.erase(pid);
    }
  }
  // send the list
  buf.write(toSendInt);
  buf.write(toSendReal);
  LOG4ESPP_INFO(theLogger, "prepared fixed pair lambda list before send particles");
}

void FixedPairLambdaList::afterRecvParticles(ParticleList &pl, InBuffer &buf) {
  std::vector<longint> receivedInt;
  std::vector<real> receivedReal;
  int n;
  longint pid1, pid2;
  real lambdaVal;
  PairsLambda::iterator it = pairsLambda.begin();
  // receive the bond list
  buf.read(receivedInt);
  buf.read(receivedReal);
  int size = receivedInt.size();
  int i = 0;
  int j = 0;
  while (i < size) {
    // unpack the list
    pid1 = receivedInt[i++];
    n = receivedInt[i++];
    for (; n > 0; --n) {
      pid2 = receivedInt[i++];
      lambdaVal = receivedReal[j++];
      it = pairsLambda.insert(it, std::make_pair(pid1, std::make_pair(pid2, lambdaVal)));
    }
  }
  if (i != size) {
    LOG4ESPP_ERROR(theLogger, "ATTETNTION:  recv particles might have read garbage.");
  }
  LOG4ESPP_INFO(theLogger, "received fixed pair lambda list after receive particles.");
}

void FixedPairLambdaList::onParticlesChanged() {
  LOG4ESPP_INFO(theLogger, "rebuild local bond list from global\n");

  System &system = storage->getSystemRef();
  esutil::Error err(system.comm);

  this->clear();
  longint lastpid1 = -1;
  Particle *p1;
  Particle *p2;
  for (PairsLambda::const_iterator it = pairsLambda.begin(); it != pairsLambda.end(); ++it) {
    if (it->first != lastpid1) {
      p1 = storage->lookupRealParticle(it->first);
      if (p1 == NULL) {
        std::stringstream msg;
        msg << "bond particle p1 " << it->first << " does not exists here";
        err.setException(msg.str());
      }
      lastpid1 = it->first;
    }
    p2 = storage->lookupLocalParticle(it->second.first);
    if (p2 == NULL) {
      std::stringstream msg;
      msg << "bond particle p2 " << it->second.first << " does not exists here";
      err.setException(msg.str());
    }
    this->push_back(ParticlePairLambda(p1, p2, it->second.second));
  }
  err.checkException();
  LOG4ESPP_INFO(theLogger, "regenerated local fixed pair list from global list");
}

void FixedPairLambdaList::registerPython() {
  using namespace espressopp::python;

  bool (FixedPairLambdaList::*pyAdd)(longint pid1, longint pid2)
  = &FixedPairLambdaList::add;

  class_<FixedPairLambdaList, shared_ptr<FixedPairLambdaList>, boost::noncopyable>
      ("FixedPairLambdaList", init<shared_ptr<storage::Storage>, real>())
      .def("add", pyAdd)
      .def("size", &FixedPairLambdaList::size)
      .def("getPairs", &FixedPairLambdaList::getPairs)
      .def("getPairsLambda", &FixedPairLambdaList::getPairsLambda)
      .def("getLambda", &FixedPairLambdaList::getLambda);
}
}  // end namespace espressopp
