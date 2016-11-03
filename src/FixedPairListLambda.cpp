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
#include "FixedPairListLambda.hpp"
#include "storage/Storage.hpp"
#include "bc/BC.hpp"
#include "esutil/Error.hpp"

namespace espressopp {

LOG4ESPP_LOGGER(FixedPairListLambda::theLogger, "FixedPairListLambda");

FixedPairListLambda::FixedPairListLambda(shared_ptr<storage::Storage> _storage, real lambda0)
    : storage(_storage), lambda0_(lambda0) {

  LOG4ESPP_INFO(theLogger, "construct FixedPairListLambda");

  sigBeforeSend = storage->beforeSendParticles.connect
      (boost::bind(&FixedPairListLambda::beforeSendParticles, this, _1, _2));
  sigAfterRecv = storage->afterRecvParticles.connect
      (boost::bind(&FixedPairListLambda::afterRecvParticles, this, _1, _2));
  sigOnParticlesChanged = storage->onParticlesChanged.connect
      (boost::bind(&FixedPairListLambda::onParticlesChanged, this));
}

bool FixedPairListLambda::add(longint pid1, longint pid2) {
  bool returnVal = true;
  if (pid1 > pid2)
    std::swap(pid1, pid2);

  System &system = storage->getSystemRef();
  esutil::Error err(system.comm);

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
    std::pair<PairsLambda::const_iterator, PairsLambda::const_iterator> equalRange = pairsLambda_.equal_range(pid1);
    if (equalRange.first != pairsLambda_.end()) {
      for (PairsLambda::const_iterator it = equalRange.first; it != equalRange.second && !found; ++it) {
        if (it->second.first == pid2)
          found = true;
      }
    }
    returnVal = !found;
    if (!found) {
      this->add(p1, p2);
      particlePairsLambda_.push_back(ParticlePairLambda(p1, p2, lambda0_));
      pairsLambda_.insert(equalRange.first, std::make_pair(pid1, std::make_pair(pid2, lambda0_)));
      onTupleAdded(pid1, pid2);
      LOG4ESPP_INFO(theLogger, "added fixed pair " << pid1 << "-" << pid2 << " to global pair list");
    }
  }
  return returnVal;
}

bool FixedPairListLambda::iadd(longint pid1, longint pid2) {
  bool returnVal = true;
  if (pid1 > pid2)
    std::swap(pid1, pid2);

  System &system = storage->getSystemRef();

  Particle *p1 = storage->lookupRealParticle(pid1);
  Particle *p2 = storage->lookupLocalParticle(pid2);

  if (!p1) {
    returnVal = false;
  } else if (!p2) {
    returnVal = false;
  }

  if (returnVal) {
    bool found = false;
    std::pair<PairsLambda::const_iterator, PairsLambda::const_iterator> equalRange = pairsLambda_.equal_range(pid1);
    if (equalRange.first != pairsLambda_.end()) {
      for (PairsLambda::const_iterator it = equalRange.first; it != equalRange.second && !found; ++it) {
        if (it->second.first == pid2)
          found = true;
      }
    }
    returnVal = !found;
    if (!found) {
      this->add(p1, p2);
      particlePairsLambda_.push_back(ParticlePairLambda(p1, p2, lambda0_));
      pairsLambda_.insert(equalRange.first, std::make_pair(pid1, std::make_pair(pid2, lambda0_)));
      onTupleAdded(pid1, pid2);
      LOG4ESPP_INFO(theLogger, "added fixed pair " << pid1 << "-" << pid2 << " to global pair list");
    }
  }
  return returnVal;
}

bool FixedPairListLambda::remove(longint pid1, longint pid2, bool no_signal) {
  LOG4ESPP_DEBUG(theLogger, "FPL remove " << pid1 << "-" << pid2);
  bool returnValue = false;
  std::pair<PairsLambda::iterator, PairsLambda::iterator> equalRange = pairsLambda_.equal_range(pid1);
  if (equalRange.first != pairsLambda_.end()) {
    for (PairsLambda::iterator it = equalRange.first; it != equalRange.second; ++it) {
      if (it->second.first == pid2) {
        LOG4ESPP_DEBUG(theLogger, "FPL, found " << it->first << " - " << it->second.first);
        if (!no_signal)
          onTupleRemoved(pid1, pid2);
        it = pairsLambda_.erase(it);
        returnValue = true;
      } else {
        it++;
      }
    }
  }
  return returnValue;
}

bool FixedPairListLambda::removeByPid1(longint pid1, bool noSignal, bool removeAll, longint removeCounter) {
  bool returnValue = false;
  std::pair<PairsLambda::iterator, PairsLambda::iterator> equalRange = pairsLambda_.equal_range(pid1);
  if (equalRange.first == pairsLambda_.end())
    return returnValue;

  if (removeAll) {
    for(PairsLambda::iterator it = equalRange.first; it != equalRange.second;) {
      if (!noSignal)
        onTupleRemoved(it->first, it->second.first);
      it = pairsLambda_.erase(it);
      returnValue = true;
    }
  } else {
    longint num_removed = 0;
    for(PairsLambda::iterator it = equalRange.first; it != equalRange.second && num_removed < removeCounter;) {
      if (!noSignal)
        onTupleRemoved(it->first, it->second.first);
      it = pairsLambda_.erase(it);
      returnValue = true;
      num_removed++;
    }
  }
  return returnValue;
}

void FixedPairListLambda::beforeSendParticles(ParticleList &pl, OutBuffer &buf) {
  std::vector<longint> toSend;
  std::vector<real> toSendLambda;
  for (ParticleList::Iterator pit(pl); pit.isValid(); ++pit) {
    longint pid = pit->id();

    LOG4ESPP_DEBUG(theLogger, "send particle with pid " << pid << ", find pairs");

    int n = globalPairs.count(pid);

    if (n > 0) {
      std::pair<PairsLambda::iterator, PairsLambda::iterator> rangeLambda = pairsLambda_.equal_range(pid);
      toSend.reserve(toSend.size() + n + 1);
      toSendLambda.reserve(toSendLambda.size() + n);
      toSend.push_back(pid);
      toSend.push_back(n);
      for (PairsLambda::const_iterator it = rangeLambda.first; it != rangeLambda.second; ++it) {
        toSend.push_back(it->second.first);
        toSendLambda.push_back(it->second.second);
        LOG4ESPP_DEBUG(theLogger, "send global bond: pid " << pid << " and partner " << it->second.first);
      }
      pairsLambda_.erase(rangeLambda.first, rangeLambda.second);
    }
  }
  buf.write(toSend);
  buf.write(toSendLambda);
  LOG4ESPP_INFO(theLogger, "prepared fixed pair list before send particles");
}

void FixedPairListLambda::afterRecvParticles(ParticleList &pl, InBuffer &buf) {
  std::vector<longint> received;
  std::vector<real> receivedLambda;
  int n;
  longint pid1, pid2;
  PairsLambda::iterator it_lambda = pairsLambda_.begin();
  buf.read(received);
  buf.read(receivedLambda);
  int size = received.size();
  int i = 0, j = 0;
  while (i < size) {
    pid1 = received[i++];
    n = received[i++];
    LOG4ESPP_DEBUG(theLogger, "recv particle " << pid1 << ", has " << n << " global pairs");
    for (; n > 0; --n) {
      pid2 = received[i++];
      real lambda = receivedLambda[j++];
      LOG4ESPP_DEBUG(theLogger, "received pair " << pid1 << " , " << pid2);
      it_lambda = pairsLambda_.insert(it_lambda, std::make_pair(pid1, std::make_pair(pid2, lambda)));
    }
  }
  if (i != size) {
    LOG4ESPP_ERROR(theLogger, "ATTETNTION:  recv particles might have read garbage\n");
  }
  LOG4ESPP_INFO(theLogger, "received fixed pair list after receive particles");
}

void FixedPairListLambda::onParticlesChanged() {
  LOG4ESPP_INFO(theLogger, "rebuild local bond list from global\n");

  System& system = storage->getSystemRef();
  esutil::Error err(system.comm);

  this->clear();
  particlePairsLambda_.clear();

  longint lastpid1 = -1;
  Particle *p1;
  Particle *p2;
  for (PairsLambda::const_iterator it = pairsLambda_.begin(); it != pairsLambda_.end(); ++it) {
    if (it->first != lastpid1) {
      p1 = storage->lookupRealParticle(it->first);
      if (p1 == NULL) {
        std::stringstream msg;
        msg << "onParticlesChanged error. Fixed Pair List particle p1 " << it->first << " does not exists here.";
        msg << " pair: " << it->first << "-" << it->second.first;
        err.setException( msg.str() );
      }
      lastpid1 = it->first;
    }
    p2 = storage->lookupLocalParticle(it->second.first);
    if (p2 == NULL) {
      std::stringstream msg;
      msg << "onParticlesChanged error. Fixed Pair List particle p2 " << it->second.first << " does not exists here.";
      msg << " p1: " << *p1;
      msg << " pair: " << it->first << "-" << it->second.first;
      err.setException( msg.str() );
    }
    this->add(p1, p2);
    particlePairsLambda_.push_back(ParticlePairLambda(p1, p2, it->second.second));
  }

  err.checkException();

  LOG4ESPP_INFO(theLogger, "regenerated local fixed pair list from global list");
}

void FixedPairListLambda::updateParticlesStorage() {
  LOG4ESPP_INFO(theLogger, "rebuild local bond list from global\n");

  System& system = storage->getSystemRef();

  this->clear();
  particlePairsLambda_.clear();
  longint lastpid1 = -1;
  Particle *p1;
  Particle *p2;
  for (PairsLambda::const_iterator it = pairsLambda_.begin(); it != pairsLambda_.end(); ++it) {
    if (it->first != lastpid1) {
      p1 = storage->lookupRealParticle(it->first);
      if (p1 == NULL) {
        std::stringstream msg;
        msg << "updateParticlesStorage error. Fixed Pair List particle p1 " << it->first << " does not exists here.";
        msg << " p1: " << *p1;
        msg << " pair: " << it->first << "-" << it->second.first;
        throw std::runtime_error(msg.str());
      }
      lastpid1 = it->first;
    }
    p2 = storage->lookupLocalParticle(it->second.first);
    if (p2 == NULL) {
      std::stringstream msg;
      msg << "updateParticlesStorage error. Fixed Pair List particle p2 " << it->second.first << " does not exists here.";
      msg << " p1: " << *p1;
      msg << " pair: " << it->first << "-" << it->second.first;
      throw std::runtime_error(msg.str());
    }
    this->add(p1, p2);
    particlePairsLambda_.push_back(ParticlePairLambda(p1, p2, it->second.second));
  }
  LOG4ESPP_INFO(theLogger, "regenerated local fixed pair list from global list");
}

real FixedPairListLambda::getLambda(longint pid1, longint pid2) {
  std::pair<PairsLambda::iterator, PairsLambda::iterator> rangeLambda = pairsLambda_.equal_range(pid1);
  for (PairsLambda::iterator it2 = rangeLambda.first; rangeLambda.second != it2; ++it2) {
    if (it2->second.first == pid2)
      return it2->second.second;
  }
}

void FixedPairListLambda::setLambda(longint pid1, longint pid2, real lambda) {
  std::pair<PairsLambda::iterator, PairsLambda::iterator> rangeLambda = pairsLambda_.equal_range(pid1);
  for (PairsLambda::iterator it2 = rangeLambda.first; rangeLambda.second != it2; ++it2) {
    if (it2->second.first == pid2)
      it2->second.second = lambda;
  }
  for (ParticlePairsLambda::iterator it = particlePairsLambda_.begin(); it != particlePairsLambda_.end(); ++it) {
    if ((it->p1->id() == pid1 && it->p2->id() == pid2) || (it->p1->id() == pid1 && it->p2->id() == pid2)) {
      it->lambda = lambda;
    }
  }
}

void FixedPairListLambda::setAllLambda(real lambda) {
  for (PairsLambda::iterator it2 = pairsLambda_.begin(); it2 != pairsLambda_.end(); ++it2) {
      it2->second.second = lambda;
  }
  for (ParticlePairsLambda::iterator it = particlePairsLambda_.begin(); it != particlePairsLambda_.end(); ++it) {
    it->lambda = lambda;
  }
}

void FixedPairListLambda::incrementAllLambda(real d_lambda) {
  for (PairsLambda::iterator it2 = pairsLambda_.begin(); it2 != pairsLambda_.end(); ++it2) {
    it2->second.second += d_lambda;
    if (it2->second.second > 1.0)
      it2->second.second = 1.0;
  }
  for (ParticlePairsLambda::iterator it = particlePairsLambda_.begin(); it != particlePairsLambda_.end(); ++it) {
    it->lambda += d_lambda;
    if (it->lambda > 1.0)
      it->lambda = 1.0;
  }
}

std::vector<longint> FixedPairListLambda::getPairList() {
  std::vector<longint> ret;
  for (PairsLambda::const_iterator it = pairsLambda_.begin(); it != pairsLambda_.end(); it++) {
    ret.push_back(it->first);
    ret.push_back(it->second.first);
  }
  return ret;
}
python::list FixedPairListLambda::getBonds() {
  python::list bonds;
  for (PairsLambda::const_iterator it=pairsLambda_.begin(); it != pairsLambda_.end(); it++) {
    bonds.append(python::make_tuple(it->first, it->second.first));
  }

  return bonds;
}
python::list FixedPairListLambda::getAllBonds() {
  std::vector<longint> local_bonds;
  std::vector<std::vector<longint> > global_bonds;
  python::list bonds;

  for (PairsLambda::const_iterator it = pairsLambda_.begin(); it != pairsLambda_.end(); it++) {
    local_bonds.push_back(it->first);
    local_bonds.push_back(it->second.first);
  }
  System& system = storage->getSystemRef();
  if (system.comm->rank() == 0) {
    mpi::gather(*system.comm, local_bonds, global_bonds, 0);
    python::tuple bond;

    for (std::vector<std::vector<longint> >::iterator it = global_bonds.begin();
         it != global_bonds.end(); ++it) {
      for (std::vector<longint>::iterator iit = it->begin(); iit != it->end();) {
        longint pid1 = *(iit++);
        longint pid2 = *(iit++);
        bonds.append(python::make_tuple(pid1, pid2));
      }
    }
  } else {
    mpi::gather(*system.comm, local_bonds, global_bonds, 0);
  }
  return bonds;
}
int FixedPairListLambda::totalSize() {
  int local_size = pairsLambda_.size();
  int global_size;
  System& system = storage->getSystemRef();
  mpi::all_reduce(*system.comm, local_size, global_size, std::plus<int>());
  return global_size;
}

void FixedPairListLambda::registerPython() {
  using namespace espressopp::python;

  bool (FixedPairListLambda::*pyAdd)(longint pid1, longint pid2) = &FixedPairListLambda::add;

  class_<FixedPairListLambda, shared_ptr<FixedPairListLambda>, boost::noncopyable >
      ("FixedPairListLambda", init<shared_ptr<storage::Storage>, real >())
      .add_property("lambda0", make_getter(&FixedPairListLambda::lambda0_), make_setter(&FixedPairListLambda::lambda0_))
      .def("add", pyAdd)
      .def("remove", &FixedPairListLambda::remove)
      .def("size", &FixedPairListLambda::size)
      .def("totalSize", &FixedPairListLambda::totalSize)
      .def("getBonds", &FixedPairListLambda::getBonds)
      .def("getAllBonds", &FixedPairListLambda::getAllBonds)
      .def("resetLongtimeMaxBondSqr", &FixedPairListLambda::resetLongtimeMaxBondSqr)
      .def("getLongtimeMaxBondSqr", &FixedPairListLambda::getLongtimeMaxBondSqr)
      .def("getLambda", &FixedPairListLambda::getLambda)
      .def("setLambda", &FixedPairListLambda::setLambda)
      .def("setAllLambda", &FixedPairListLambda::setAllLambda);
}
}  // end namespace espressopp
