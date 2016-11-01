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

#ifndef _FIXEDTRIPLELAMBDALIST_HPP
#define _FIXEDTRIPLELAMBDALIST_HPP

#include "log4espp.hpp"
#include "python.hpp"
#include "FixedTripleList.hpp"
#include "Particle.hpp"
#include "esutil/ESPPIterator.hpp"
#include <map>
#include <boost/signals2.hpp>
#include <boost/unordered_map.hpp>
#include "types.hpp"

namespace espressopp {

class ParticleTripleLambda {
 public:
  ParticleTripleLambda(Particle *p1, Particle *p2, Particle *p3, real l)
      : first(p1), second(p2), third(p3), lambda(l) {}

  ParticleTripleLambda(Particle &p1, Particle &p2, Particle &p3, real l)
      : first(&p1), second(&p2), third(&p3), lambda(l) {}

  Particle *first;
  Particle *second;
  Particle *third;
  real lambda;
};

class FixedTripleListLambda: public FixedTripleList {
 protected:
  typedef boost::unordered_multimap<longint, std::pair <longint, std::pair<longint, real> > > TriplesLambda;
  boost::signals2::connection sigAfterRecv, sigOnParticleChanged, sigBeforeSend;
  shared_ptr <storage::Storage> storage;
  using TripleList::add;

 public:
  typedef std::vector<ParticleTripleLambda> ParticleTriplesLambda;
  typedef esutil::ESPPIterator<std::vector<ParticleTripleLambda> > IteratorParticleLambda;
  FixedTripleListLambda(shared_ptr <storage::Storage> _storage, real initLambda);

  bool add(longint pid1, longint pid2, longint pid3);
  bool iadd(longint pid1, longint pid2, longint pid3);
  void beforeSendParticles(ParticleList &pl, class OutBuffer &buf);
  void afterRecvParticles(ParticleList &pl, class InBuffer &buf);
  void onParticlesChanged();
  bool remove(longint pid1, longint pid2, longint pid3);
  void updateParticlesStorage();

  python::list getTriples();
  python::list getTriplesLambda();
  /** Get the number of triples in the GlobalTriples list */
  virtual int size() { return triplesLambda_.size(); }
  virtual int totalSize();

  real getLambda(longint pid1, longint pid2, longint pid3);
  void setLambda(longint pid1, longint pid2, longint pid3, real lambda);
  void setAllLambda(real lambda);
  void incrementAllLambda(real d_lambda);

  ParticleTriplesLambda& getParticleTriples() { return particleTriplesLambda_; }

  boost::signals2::signal3 <void, longint, longint, longint> onTupleAdded;
  boost::signals2::signal3 <void, longint, longint, longint> onTupleRemoved;

  static void registerPython();
 private:
  real lambda0_;
  TriplesLambda triplesLambda_;
  ParticleTriplesLambda particleTriplesLambda_;
  static LOG4ESPP_DECL_LOGGER(theLogger);

};
}

#endif

