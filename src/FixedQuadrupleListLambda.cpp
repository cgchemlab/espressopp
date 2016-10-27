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
#include "FixedQuadrupleListLambda.hpp"
#include "storage/Storage.hpp"

#include "esutil/Error.hpp"

namespace espressopp {

LOG4ESPP_LOGGER(FixedQuadrupleListLambda::theLogger, "FixedQuadrupleListLambda");

FixedQuadrupleListLambda::FixedQuadrupleListLambda(shared_ptr<storage::Storage> _storage, real initLambda)
    : storage(_storage), globalQuadruples(), initLambda_(initLambda) {
  LOG4ESPP_INFO(theLogger, "construct FixedQuadrupleListLambda");

  sigBeforeSend = storage->beforeSendParticles.connect
      (boost::bind(&FixedQuadrupleListLambda::beforeSendParticles, this, _1, _2));
  sigAfterRecv = storage->afterRecvParticles.connect
      (boost::bind(&FixedQuadrupleListLambda::afterRecvParticles, this, _1, _2));
  sigOnParticlesChanged = storage->onParticlesChanged.connect
      (boost::bind(&FixedQuadrupleListLambda::onParticlesChanged, this));
}

FixedQuadrupleListLambda::~FixedQuadrupleListLambda() {
  LOG4ESPP_INFO(theLogger, "~FixedQuadrupleListLambda");

  sigBeforeSend.disconnect();
  sigAfterRecv.disconnect();
  sigOnParticlesChanged.disconnect();
}

bool FixedQuadrupleListLambda::iadd(longint pid1, longint pid2, longint pid3, longint pid4) {
  // here we assume pid1 < pid2 < pid3 < pid4
  bool returnVal = true;
  System &system = storage->getSystemRef();

  // ADD THE LOCAL QUADRUPLET
  Particle *p1 = storage->lookupRealParticle(pid1);
  Particle *p2 = storage->lookupLocalParticle(pid2);
  Particle *p3 = storage->lookupLocalParticle(pid3);
  Particle *p4 = storage->lookupLocalParticle(pid4);
  if (!p1) {
    // Particle does not exist here, return false
    returnVal = false;
  } else {
    std::stringstream msg;
    if (!p2) {
      msg << "quadruple particle p2 " << pid2 << " does not exists here and cannot be added";
      throw std::runtime_error(msg.str());
    }
    if (!p3) {
      msg << "quadruple particle p3 " << pid3 << " does not exists here and cannot be added";
      throw std::runtime_error(msg.str());
    }
    if (!p4) {
      msg << "quadruple particle p4 " << pid4 << " does not exists here and cannot be added";
      throw std::runtime_error(msg.str());
    }
  }

  if (returnVal) {
    bool found = false;
    std::pair<GlobalQuadruples::const_iterator,
              GlobalQuadruples::const_iterator> equalRange
        = globalQuadruples.equal_range(pid1);
    if (equalRange.first != globalQuadruples.end()) {
      for (GlobalQuadruples::const_iterator it = equalRange.first; it != equalRange.second && !found; ++it)
        if (it->second.first == Triple<longint, longint, longint>(pid2, pid3, pid4))
          found = true;
    }
    returnVal = !found;
    if (!found) {
      // add the quadruple locally
      this->push_back(ParticleQuadrupleLambda(p1, p2, p3, p4, initLambda_));
      // if not, insert the new quadruple
      globalQuadruples.insert(equalRange.first,
                              std::make_pair(pid1,
                                             std::make_pair(
                                                 Triple<longint, longint, longint>(pid2, pid3, pid4),
                                                 initLambda_)));
      onTupleAdded(pid1, pid2, pid3, pid4);
      LOG4ESPP_INFO(theLogger, "added fixed quadruple to global quadruple list");
    }
  }
  return returnVal;
}


bool FixedQuadrupleListLambda::add(longint pid1, longint pid2, longint pid3, longint pid4) {
  // here we assume pid1 < pid2 < pid3 < pid4
  bool returnVal = true;
  System &system = storage->getSystemRef();
  esutil::Error err(system.comm);

  // ADD THE LOCAL QUADRUPLET
  Particle *p1 = storage->lookupRealParticle(pid1);
  Particle *p2 = storage->lookupLocalParticle(pid2);
  Particle *p3 = storage->lookupLocalParticle(pid3);
  Particle *p4 = storage->lookupLocalParticle(pid4);
  if (!p1) {
    // Particle does not exist here, return false
    returnVal = false;
  } else {
    if (!p2) {
      std::stringstream msg;
      msg << "quadruple particle p2 " << pid2 << " does not exists here and cannot be added";
      err.setException(msg.str());
    }
    if (!p3) {
      std::stringstream msg;
      msg << "quadruple particle p3 " << pid3 << " does not exists here and cannot be added";
      err.setException(msg.str());
    }
    if (!p4) {
      std::stringstream msg;
      msg << "quadruple particle p4 " << pid4 << " does not exists here and cannot be added";
      err.setException(msg.str());
    }
  }
  err.checkException();

  if (returnVal) {
    bool found = false;
    std::pair<GlobalQuadruples::const_iterator,
              GlobalQuadruples::const_iterator> equalRange
        = globalQuadruples.equal_range(pid1);
    if (equalRange.first != globalQuadruples.end()) {
      for (GlobalQuadruples::const_iterator it = equalRange.first;
           it != equalRange.second && !found; ++it)
        if (it->second.first == Triple<longint, longint, longint>(pid2, pid3, pid4))
          found = true;
    }
    returnVal = !found;
    if (!found) {
      this->push_back(ParticleQuadrupleLambda(p1, p2, p3, p4, initLambda_));
      globalQuadruples.insert(equalRange.first,
                              std::make_pair(pid1,
                                             std::make_pair(
                                                 Triple<longint, longint, longint>(pid2, pid3, pid4),
                                                 initLambda_)));
      onTupleAdded(pid1, pid2, pid3, pid4);
      LOG4ESPP_INFO(theLogger, "added fixed quadruple to global quadruple list");
    }
  }
  return returnVal;
}

bool FixedQuadrupleListLambda::remove(longint pid1, longint pid2, longint pid3, longint pid4) {
  // here we assume pid1 < pid2 < pid3 < pid4
  bool found = true;
  System &system = storage->getSystemRef();

  // ADD THE LOCAL QUADRUPLET
  Particle *p1 = storage->lookupRealParticle(pid1);
  Particle *p2 = storage->lookupLocalParticle(pid2);
  Particle *p3 = storage->lookupLocalParticle(pid3);
  Particle *p4 = storage->lookupLocalParticle(pid4);
  if (!p1) {
    // Particle does not exist here, return false
    found = false;
  } else {
    std::stringstream msg;
    if (!p2) {
      msg << "quadruple particle p2 " << pid2 << " does not exists here and cannot be added";
      throw std::runtime_error(msg.str());
    }
    if (!p3) {
      msg << "quadruple particle p3 " << pid3 << " does not exists here and cannot be added";
      throw std::runtime_error(msg.str());
    }
    if (!p4) {
      msg << "quadruple particle p4 " << pid4 << " does not exists here and cannot be added";
      throw std::runtime_error(msg.str());
    }
  }

  bool returnVal = false;
  if (found) {
    std::pair<GlobalQuadruples::iterator, GlobalQuadruples::iterator> equalRange
        = globalQuadruples.equal_range(pid1);
    if (equalRange.first != globalQuadruples.end()) {
      // otherwise test whether the quadruple already exists
      for (GlobalQuadruples::iterator it = equalRange.first; it != equalRange.second;)
        if (it->second.first == Triple<longint, longint, longint>(pid2, pid3, pid4)) {
          onTupleRemoved(pid1, pid2, pid3, pid4);
          returnVal = true;
          it = globalQuadruples.erase(it);
        } else {
          ++it;
        }
    }
  }
  return returnVal;
}

python::list FixedQuadrupleListLambda::getQuadruples() {
  python::tuple quadruple;
  python::list quadruples;
  for (GlobalQuadruples::const_iterator it = globalQuadruples.begin(); it != globalQuadruples.end(); it++) {
    quadruple = python::make_tuple(it->first, it->second.first.first, it->second.first.second, it->second.first.third);
    quadruples.append(quadruple);
  }

  return quadruples;
}

void FixedQuadrupleListLambda::beforeSendParticles(ParticleList &pl, OutBuffer &buf) {
  std::vector<longint> toSend;
  std::vector<real> toSendReal;
  for (ParticleList::Iterator pit(pl); pit.isValid(); ++pit) {
    longint pid = pit->id();
    int n = globalQuadruples.count(pid);

    if (n > 0) {
      std::pair<GlobalQuadruples::const_iterator, GlobalQuadruples::const_iterator>
          equalRange = globalQuadruples.equal_range(pid);
      toSend.reserve(toSend.size() + 3 * n + 1);
      toSend.push_back(pid);
      toSend.push_back(n);
      for (GlobalQuadruples::const_iterator it = equalRange.first; it != equalRange.second; ++it) {
        toSend.push_back(it->second.first.first);
        toSend.push_back(it->second.first.second);
        toSend.push_back(it->second.first.third);
        toSendReal.push_back(it->second.second);
      }
      // delete all of these quadruples from the global list
      globalQuadruples.erase(equalRange.first, equalRange.second);
    }
  }
  // send the list
  buf.write(toSend);
  buf.write(toSendReal);
  LOG4ESPP_INFO(theLogger, "prepared fixed quadruple list before send particles");
}

void FixedQuadrupleListLambda::afterRecvParticles(ParticleList &pl, InBuffer &buf) {

  std::vector<longint> received;
  std::vector<real> receivedReal;
  int n;
  real lambdaVal;
  longint pid1, pid2, pid3, pid4;
  GlobalQuadruples::iterator it = globalQuadruples.begin();
  buf.read(received);
  buf.read(receivedReal);
  int size = received.size();
  int i = 0, j = 0;
  while (i < size) {
    pid1 = received[i++];
    n = received[i++];
    for (; n > 0; --n) {
      pid2 = received[i++];
      pid3 = received[i++];
      pid4 = received[i++];
      lambdaVal = receivedReal[j++];
      it = globalQuadruples.insert(it, std::make_pair(pid1,
                                                      std::make_pair(
                                                          Triple<longint, longint, longint>(pid2, pid3, pid4),
                                                          lambdaVal)));
    }
  }
  if (i != size) {
    printf("ATTETNTION:  recv particles might have read garbage\n");
  }
  LOG4ESPP_INFO(theLogger, "received fixed quadruple list after receive particles");
}

void FixedQuadrupleListLambda::onParticlesChanged() {
  System &system = storage->getSystemRef();
  esutil::Error err(system.comm);

  this->clear();
  longint lastpid1 = -1;
  Particle *p1;
  Particle *p2;
  Particle *p3;
  Particle *p4;
  for (GlobalQuadruples::const_iterator it = globalQuadruples.begin(); it != globalQuadruples.end(); ++it) {
    if (it->first != lastpid1) {
      p1 = storage->lookupRealParticle(it->first);
      if (p1 == NULL) {
        std::stringstream msg;
        msg << "quadruple particle p1 " << it->first << " does not exists here";
        err.setException(msg.str());
      }
      lastpid1 = it->first;
    }
    p2 = storage->lookupLocalParticle(it->second.first.first);
    if (p2 == NULL) {
      std::stringstream msg;
      msg << "quadruple particle p2 " << it->second.first.first << " does not exists here";
      err.setException(msg.str());
    }
    p3 = storage->lookupLocalParticle(it->second.first.second);
    if (p3 == NULL) {
      std::stringstream msg;
      msg << "quadruple particle p3 " << it->second.first.second << " does not exists here";
      err.setException(msg.str());
    }
    p4 = storage->lookupLocalParticle(it->second.first.third);
    if (p4 == NULL) {
      std::stringstream msg;
      msg << "quadruple particle p4 " << it->second.first.third << " does not exists here";
      err.setException(msg.str());
    }
    this->push_back(ParticleQuadrupleLambda(p1, p2, p3, p4, it->second.second));
  }
  LOG4ESPP_INFO(theLogger, "regenerated local fixed quadruple list from global list");
}

void FixedQuadrupleListLambda::updateParticlesStorage() {
  System &system = storage->getSystemRef();

  this->clear();
  longint lastpid1 = -1;
  Particle *p1;
  Particle *p2;
  Particle *p3;
  Particle *p4;
  for (GlobalQuadruples::const_iterator it = globalQuadruples.begin(); it != globalQuadruples.end(); ++it) {
    if (it->first != lastpid1) {
      p1 = storage->lookupRealParticle(it->first);
      if (p1 == NULL) {
        std::stringstream msg;
        msg << "quadruple particle p1 " << it->first << " does not exists here";
        throw std::runtime_error(msg.str());
      }
      lastpid1 = it->first;
    }
    p2 = storage->lookupLocalParticle(it->second.first.first);
    if (p2 == NULL) {
      std::stringstream msg;
      msg << "quadruple particle p2 " << it->second.first.first << " does not exists here";
      throw std::runtime_error(msg.str());
    }
    p3 = storage->lookupLocalParticle(it->second.first.second);
    if (p3 == NULL) {
      std::stringstream msg;
      msg << "quadruple particle p3 " << it->second.first.second << " does not exists here";
      throw std::runtime_error(msg.str());
    }
    p4 = storage->lookupLocalParticle(it->second.first.third);
    if (p4 == NULL) {
      std::stringstream msg;
      msg << "quadruple particle p4 " << it->second.first.third << " does not exists here";
      throw std::runtime_error(msg.str());
    }
    this->push_back(ParticleQuadrupleLambda(p1, p2, p3, p4, it->second.second));
  }
  LOG4ESPP_INFO(theLogger, "regenerated local fixed quadruple list from global list");
}

int FixedQuadrupleListLambda::totalSize() {
  int local_size = globalQuadruples.size();
  int global_size;
  System &system = storage->getSystemRef();
  mpi::all_reduce(*system.comm, local_size, global_size, std::plus<int>());
  return global_size;
}

python::list FixedQuadrupleListLambda::getQuadruplesLambda() {
  python::list quaruples;
  for (GlobalQuadruples::const_iterator it=globalQuadruples.begin(); it != globalQuadruples.end(); it++) {
    quaruples.append(python::make_tuple(
        it->first,
        it->second.first.first,
        it->second.first.second,
        it->second.first.third));
  }
  return quaruples;
}

real FixedQuadrupleListLambda::getLambda(longint pid1, longint pid2, longint pid3, longint pid4) {
  real returnVal = -3;

  bool found = false;
  std::pair<GlobalQuadruples::const_iterator,
            GlobalQuadruples::const_iterator> equalRange = globalQuadruples.equal_range(pid1);
  if (equalRange.first != globalQuadruples.end()) {
    for (GlobalQuadruples::const_iterator it = equalRange.first; it != equalRange.second; ++it) {
      if (it->second.first == Triple<longint, longint, longint>(pid2, pid3, pid4)) {
        return it->second.second;
      }
    }
  }

  return returnVal;
}

void FixedQuadrupleListLambda::setLambda(longint pid1, longint pid2, longint pid3, longint pid4, real lambda) {
  std::pair<GlobalQuadruples::iterator,
            GlobalQuadruples::iterator> equalRange = globalQuadruples.equal_range(pid1);
  if (equalRange.first != globalQuadruples.end()) {
    for (GlobalQuadruples::iterator it = equalRange.first; it != equalRange.second; ++it) {
      if (it->second.first == Triple<longint, longint, longint>(pid2, pid3, pid4)) {
        it->second.second = lambda;
      }
    }
  }

  for (FixedQuadrupleListLambda::Iterator it(*this); it.isValid(); ++it) {
    const Particle &p1 = *it->first;
    const Particle &p2 = *it->second;
    const Particle &p3 = *it->third;
    const Particle &p4 = *it->fourth;
    if ((p1.id() == pid1 && p2.id() == pid2 && p3.id() == pid3 && p4.id() == pid4) ||
        (p1.id() == pid4 && p2.id() == pid3 && p3.id() == pid2 && p4.id() == pid1) ){
      it->fifth = lambda;
    }
  }
}

void FixedQuadrupleListLambda::setLambdaAll(real lambda) {
  for (GlobalQuadruples::iterator it = globalQuadruples.begin(); it != globalQuadruples.end(); ++it) {
    it->second.second = lambda;
  }

  for (FixedQuadrupleListLambda::Iterator it(*this); it.isValid(); ++it) {
    it->fifth = lambda;
  }
}

void FixedQuadrupleListLambda::registerPython() {
  using namespace espressopp::python;

  bool (FixedQuadrupleListLambda::*pyAdd)(longint pid1, longint pid2,
                                          longint pid3, longint pid4) = &FixedQuadrupleListLambda::add;

  class_<FixedQuadrupleListLambda, shared_ptr<FixedQuadrupleListLambda>, boost::noncopyable>
      ("FixedQuadrupleListLambda", init<shared_ptr<storage::Storage>, real >())
      .def("add", pyAdd)
      .def("size", &FixedQuadrupleListLambda::size)
      .def("totalSize", &FixedQuadrupleListLambda::totalSize)
      .def("getQuadruples", &FixedQuadrupleListLambda::getQuadruples)
      .def("getQuadruplesLambda", &FixedQuadrupleListLambda::getQuadruplesLambda)
      .def("getLambda", &FixedQuadrupleListLambda::getLambda)
      .def("setLambda", &FixedQuadrupleListLambda::setLambda)
      .def("setLambdaAll", &FixedQuadrupleListLambda::setLambdaAll);
}
}  // end namespace espressopp
