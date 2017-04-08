/*
  Copyright (C) 2016-2017
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
#include "FixedTripleListLambda.hpp"
#include "storage/Storage.hpp"
#include "bc/BC.hpp"
#include "esutil/Error.hpp"

namespace espressopp {

LOG4ESPP_LOGGER(FixedTripleListLambda::theLogger, "FixedTripleListLambda");

FixedTripleListLambda::FixedTripleListLambda(shared_ptr<storage::Storage> _storage, real initLambda)
    : storage(_storage), lambda0_(initLambda) {
  LOG4ESPP_INFO(theLogger, "construct FixedTripleListLambda");

  sigBeforeSend = storage->beforeSendParticles.connect
      (boost::bind(&FixedTripleListLambda::beforeSendParticles, this, _1, _2));
  sigAfterRecv = storage->afterRecvParticles.connect
      (boost::bind(&FixedTripleListLambda::afterRecvParticles, this, _1, _2));
  sigOnParticleChanged = storage->onParticlesChanged.connect
      (boost::bind(&FixedTripleListLambda::onParticlesChanged, this));
}

bool FixedTripleListLambda::iadd(longint pid1, longint pid2, longint pid3) {
  bool returnVal = true;
  System& system = storage->getSystemRef();

  Particle *p1 = storage->lookupLocalParticle(pid1);
  Particle *p2 = storage->lookupRealParticle(pid2);
  Particle *p3 = storage->lookupLocalParticle(pid3);

  if (!p2){
    returnVal = false;
  } else {
    std::stringstream msg;
    if (!p1) {
      msg << "adding error: triple particle p1 " << pid1 <<
          " does not exists here and cannot be added";
      msg << " triplet: " << pid1 << "-" << pid2 << "-" << pid3;
      throw std::runtime_error(msg.str());
    }
    if (!p3) {
      msg << "adding error: triple particle p3 " << pid3 <<
          " does not exists here and cannot be added";
      msg << " triplet: " << pid1 << "-" << pid2 << "-" << pid3;
      throw std::runtime_error(msg.str());
    }
  }

  if (returnVal) {
    bool found = false;
    std::pair<TriplesLambda::const_iterator,
              TriplesLambda::const_iterator> equalRange = triplesLambda_.equal_range(pid2);
    if (equalRange.first != triplesLambda_.end()) {
      for (TriplesLambda::const_iterator it = equalRange.first;
           it != equalRange.second && !found; ++it) {
        found = ((it->second.first == pid1 && it->second.second.first == pid3) ||
            (it->second.first == pid3 && it->second.second.first == pid1));
      }
    }
    returnVal = !found;
    if (!found) {
      this->add(p1, p2, p3);
      particleTriplesLambda_.push_back(ParticleTripleLambda(p1, p2, p3, lambda0_));
      triplesLambda_.insert(equalRange.first,
                            std::make_pair(pid2,
                                           std::make_pair(pid1,
                                                          std::make_pair(pid3, lambda0_))));
      onTupleAdded(pid1, pid2, pid3);
    }
    LOG4ESPP_INFO(theLogger, "added fixed triple to global triple list");
  }
  return returnVal;
}

bool FixedTripleListLambda::add(longint pid1, longint pid2, longint pid3) {
  bool returnVal = true;
  System& system = storage->getSystemRef();
  esutil::Error err(system.comm);

  Particle *p1 = storage->lookupLocalParticle(pid1);
  Particle *p2 = storage->lookupRealParticle(pid2);
  Particle *p3 = storage->lookupLocalParticle(pid3);

  if (!p2){
    returnVal = false;
  } else {
    if (!p1) {
      std::stringstream msg;
      msg << "adding error: triple particle p1 " << pid1 << " does not exists here and cannot be added";
      msg << " triplet: " << pid1 << "-" << pid2 << "-" << pid3;
      err.setException( msg.str() );
    }
    if (!p3) {
      std::stringstream msg;
      msg << "adding error: triple particle p3 " << pid3 << " does not exists here and cannot be added";
      msg << " triplet: " << pid1 << "-" << pid2 << "-" << pid3;
      err.setException( msg.str() );
    }
  }

  err.checkException();

  if (returnVal) {
    bool found = false;
    std::pair<TriplesLambda::const_iterator,
              TriplesLambda::const_iterator> equalRange = triplesLambda_.equal_range(pid2);
    if (equalRange.first != triplesLambda_.end()) {
      for (TriplesLambda::const_iterator it = equalRange.first; it != equalRange.second && !found; ++it) {
        found = ((it->second.first == pid1 && it->second.second.first == pid3) ||
            (it->second.first == pid3 && it->second.second.first == pid1));
      }
    }
    returnVal = !found;
    if (!found) {
      this->add(p1, p2, p3);
      particleTriplesLambda_.push_back(ParticleTripleLambda(p1, p2, p3, lambda0_));
      triplesLambda_.insert(equalRange.first,
                           std::make_pair(pid2,
                                          std::make_pair(pid1,
                                                         std::make_pair(pid3, lambda0_))));
      onTupleAdded(pid1, pid2, pid3);
    }
    LOG4ESPP_INFO(theLogger, "added fixed triple to global triple list");
  }
  return returnVal;
}

bool FixedTripleListLambda::remove(longint pid1, longint pid2, longint pid3) {

  bool returnVal = false;
  // Remove entries.
  std::pair<TriplesLambda::iterator, TriplesLambda::iterator> equalRange =
      triplesLambda_.equal_range(pid2);
  if (equalRange.first != triplesLambda_.end()) {
    // otherwise test whether the triple already exists
    for (TriplesLambda::iterator it = equalRange.first; it != equalRange.second;) {
        if ((it->second.first == pid1 && it->second.second.first == pid3) ||
            (it->second.first == pid3 && it->second.second.first == pid1)) {
        LOG4ESPP_DEBUG(theLogger, "removed triple " << it->first << "-" << it->second.first
                                                    << "-" << it->second.second.first);
        onTupleRemoved(it->second.first, it->first, it->second.second.first);
        it = triplesLambda_.erase(it);
        returnVal = true;
      } else {
        ++it;
      }
    }
  }
  return returnVal;
}


bool FixedTripleListLambda::removeByBond(longint pid1, longint pid2) {
  bool return_val = false;

  for (TriplesLambda::iterator it = triplesLambda_.begin(); it != triplesLambda_.end();) {
    longint a1 = it->second.first;
    longint a2 = it->first;
    longint a3 = it->second.second.first;
    if ((a1 == pid1 && a2 == pid2) || (a2 == pid1 && a3 == pid2) ||
        (a1 == pid2 && a2 == pid1) || (a2 == pid2 && a3 == pid1)) {
      it = triplesLambda_.erase(it);
      return_val = true;
    } else {
      ++it;
    }
  }

  return return_val;
}

int FixedTripleListLambda::totalSize() {
  int local_size = triplesLambda_.size();
  int global_size;
  System& system = storage->getSystemRef();
  mpi::all_reduce(*system.comm, local_size, global_size, std::plus<int>());
  return global_size;
}

python::list FixedTripleListLambda::getTriples() {
  python::tuple triple;
  python::list triples;
  for (TriplesLambda::const_iterator it=triplesLambda_.begin(); it != triplesLambda_.end(); it++) {
    triple = python::make_tuple(it->second.first, it->first, it->second.second.first);
    triples.append(triple);
  }

  return triples;
}

python::list FixedTripleListLambda::getTriplesLambda() {
  python::tuple triple;
  python::list triples;
  for (TriplesLambda::const_iterator it=triplesLambda_.begin(); it != triplesLambda_.end(); it++) {
    triple = python::make_tuple(it->second.first, it->first, it->second.second.first, it->second.second.second);
    triples.append(triple);
  }

  return triples;
}

real FixedTripleListLambda::getLambda(longint pid1, longint pid2, longint pid3) {
  real returnVal = -3;

  bool found = false;
  std::pair<TriplesLambda::const_iterator,
            TriplesLambda::const_iterator> equalRange = triplesLambda_.equal_range(pid2);
  if (equalRange.first != triplesLambda_.end()) {
    for (TriplesLambda::const_iterator it = equalRange.first; it != equalRange.second && !found; ++it) {
      if ((it->second.first == pid1 && it->second.second.first == pid3) ||
          (it->second.first == pid3 && it->second.second.first == pid1)) {
        return it->second.second.second;
      }
    }
  }

  return returnVal;
}

void FixedTripleListLambda::setLambda(longint pid1, longint pid2, longint pid3, real lambda) {

  std::pair<TriplesLambda::iterator,
            TriplesLambda::iterator> equalRange = triplesLambda_.equal_range(pid2);
  if (equalRange.first != triplesLambda_.end()) {
    for (TriplesLambda::iterator it = equalRange.first; it != equalRange.second; ++it) {
      if ((it->second.first == pid1 && it->second.second.first == pid3) ||
          (it->second.first == pid3 && it->second.second.first == pid1)) {
        it->second.second.second = lambda;
      }
    }
  }

  for (ParticleTriplesLambda::iterator it = particleTriplesLambda_.begin(); it != particleTriplesLambda_.end(); ++it) {
    const Particle &p1 = *it->first;
    const Particle &p2 = *it->second;
    const Particle &p3 = *it->third;
    if ((p1.id() == pid1 && p2.id() == pid2 && p3.id() == pid3) ||
        (p1.id() == pid3 && p2.id() == pid2 && p3.id() == pid1) ){
      it->lambda = lambda;
    }
  }
}

void FixedTripleListLambda::setAllLambda(real lambda) {
  for (TriplesLambda::iterator it = triplesLambda_.begin(); it != triplesLambda_.end(); ++it) {
    it->second.second.second = lambda;
  }

  for (ParticleTriplesLambda::iterator it = particleTriplesLambda_.begin(); it != particleTriplesLambda_.end(); ++it) {
    it->lambda = lambda;
  }
}

void FixedTripleListLambda::incrementAllLambda(real d_lambda) {
  for (TriplesLambda::iterator it = triplesLambda_.begin(); it != triplesLambda_.end(); ++it) {
    it->second.second.second += d_lambda;
    if (it->second.second.second > 1.0)
      it->second.second.second = 1.0;
  }

  for (ParticleTriplesLambda::iterator it = particleTriplesLambda_.begin(); it != particleTriplesLambda_.end(); ++it) {
    it->lambda += d_lambda;
    if (it->lambda > 1.0)
      it->lambda = 1.0;
  }
}

void FixedTripleListLambda::beforeSendParticles(ParticleList& pl, OutBuffer& buf) {
  std::vector<longint> toSend;
  std::vector<real> toSendReal;
  // loop over the particle list
  for (ParticleList::Iterator pit(pl); pit.isValid(); ++pit) {
    longint pid = pit->id();
    int n = triplesLambda_.count(pid);
    if (n > 0) {
      std::pair<TriplesLambda::const_iterator, TriplesLambda::const_iterator> equalRange = triplesLambda_.equal_range(pid);

      toSend.reserve(toSend.size() + 2*n + 1);
      toSendReal.reserve(toSendReal.size() + n);
      toSend.push_back(pid);
      toSend.push_back(n);
      for (TriplesLambda::const_iterator it = equalRange.first;
           it != equalRange.second; ++it) {
        toSend.push_back(it->second.first);
        toSend.push_back(it->second.second.first);
        toSendReal.push_back(it->second.second.second);
      }
      triplesLambda_.erase(equalRange.first, equalRange.second);
    }
  }
  buf.write(toSend);
  buf.write(toSendReal);
  LOG4ESPP_INFO(theLogger, "prepared fixed triple list before send particles");
}

void FixedTripleListLambda::afterRecvParticles(ParticleList &pl, InBuffer& buf) {
  std::vector< longint > receivedInt;
  std::vector<real> receivedReal;
  int n;
  longint pid1, pid2, pid3;
  real lambdaVal;
  TriplesLambda::iterator it = triplesLambda_.begin();
  // receive the triple list
  buf.read(receivedInt);
  buf.read(receivedReal);
  int size = receivedInt.size();
  int i = 0, j = 0;
  while (i < size) {
    // unpack the list
    pid2 = receivedInt[i++];
    n = receivedInt[i++];
    for (; n > 0; --n) {
      pid1 = receivedInt[i++];
      pid3 = receivedInt[i++];
      lambdaVal = receivedReal[j++];
      it = triplesLambda_.insert(
          it, std::make_pair(pid2, std::make_pair(pid1, std::make_pair(pid3, lambdaVal))));
    }
  }
  if (i != size) {
    printf("ATTETNTION:  recv particles might have read garbage\n");
  }
  LOG4ESPP_INFO(theLogger, "received fixed triple list after receive particles");
}

void FixedTripleListLambda::onParticlesChanged() {

  System& system = storage->getSystemRef();
  esutil::Error err(system.comm);

  // (re-)generate the local triple list from the global list
  this->clear();
  particleTriplesLambda_.clear();
  longint lastpid2 = -1;
  Particle *p1;
  Particle *p2;
  Particle *p3;
  for (TriplesLambda::const_iterator it = triplesLambda_.begin();
       it != triplesLambda_.end(); ++it) {
    //printf("lookup global triple %d %d %d\n", it->first, it->second.first, it->second.second);
    if (it->first != lastpid2) {
      p2 = storage->lookupRealParticle(it->first);
      if (p2 == NULL) {
        std::stringstream msg;
        msg << "triple particle p2 " << it->first << " does not exists here";
        err.setException( msg.str() );
      }
      lastpid2 = it->first;
    }
    p1 = storage->lookupLocalParticle(it->second.first);
    if (p1 == NULL) {
      std::stringstream msg;
      msg << "triple particle p1 " << it->second.first << " does not exists here";
      err.setException( msg.str() );
    }
    p3 = storage->lookupLocalParticle(it->second.second.first);
    if (p3 == NULL) {
      std::stringstream msg;
      msg << "triple particle p3 " << it->second.second.first << " does not exists here";
      err.setException( msg.str() );
    }
    this->add(p1, p2, p3);
    particleTriplesLambda_.push_back(ParticleTripleLambda(p1, p2, p3, it->second.second.second));
  }
  err.checkException();

  LOG4ESPP_INFO(theLogger, "regenerated local fixed triple list from global list");
}

void FixedTripleListLambda::updateParticlesStorage() {
  System& system = storage->getSystemRef();

  this->clear();
  particleTriplesLambda_.clear();
  longint lastpid2 = -1;
  Particle *p1;
  Particle *p2;
  Particle *p3;
  for (TriplesLambda::const_iterator it = triplesLambda_.begin(); it != triplesLambda_.end(); ++it) {
    if (it->first != lastpid2) {
      p2 = storage->lookupRealParticle(it->first);
      if (p2 == NULL) {
        std::stringstream msg;
        msg << "triple particle p2 " << it->first << " does not exists here";
        throw std::runtime_error(msg.str());
      }
      lastpid2 = it->first;
    }
    p1 = storage->lookupLocalParticle(it->second.first);
    if (p1 == NULL) {
      std::stringstream msg;
      msg << "triple particle p1 " << it->second.first << " does not exists here";
      throw std::runtime_error(msg.str());
    }
    p3 = storage->lookupLocalParticle(it->second.second.first);
    if (p3 == NULL) {
      std::stringstream msg;
      msg << "triple particle p3 " << it->second.second.first << " does not exists here";
      throw std::runtime_error(msg.str());
    }
    this->add(p1, p2, p3);
    particleTriplesLambda_.push_back(ParticleTripleLambda(p1, p2, p3, it->second.second.second));
  }

  LOG4ESPP_INFO(theLogger, "regenerated local fixed triple list from global list");
}

void FixedTripleListLambda::registerPython() {
  using namespace espressopp::python;

  bool (FixedTripleListLambda::*pyAdd)(longint pid1, longint pid2, longint pid3) = &FixedTripleListLambda::add;

  class_<FixedTripleListLambda, shared_ptr<FixedTripleListLambda>, boost::noncopyable>
      ("FixedTripleListLambda", init<shared_ptr<storage::Storage>, real>())
      .add_property("lambda0", make_getter(&FixedTripleListLambda::lambda0_),
                    make_setter(&FixedTripleListLambda::lambda0_))
      .def("add", pyAdd)
      .def("size", &FixedTripleListLambda::size)
      .def("getTriples", &FixedTripleListLambda::getTriples)
      .def("getTriplesLambda", &FixedTripleListLambda::getTriplesLambda)
      .def("getLambda", &FixedTripleListLambda::getLambda)
      .def("setLambda", &FixedTripleListLambda::setLambda)
      .def("setAllLambda", &FixedTripleListLambda::setAllLambda);
}


}  // end namespace espressopp
