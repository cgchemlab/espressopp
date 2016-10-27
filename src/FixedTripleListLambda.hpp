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
#include "Particle.hpp"
#include "esutil/ESPPIterator.hpp"
#include <map>
#include <boost/signals2.hpp>
#include <boost/unordered_map.hpp>
#include "types.hpp"

namespace espressopp {

class ParticleTripleLambda : public Quadruple<class Particle*, class Particle*, class Particle*, real> {
 private:
  typedef Quadruple<class Particle*, class Particle*, class Particle*, real> Super;
 public:
  ParticleTripleLambda() : Super() {}
  ParticleTripleLambda(Particle* p1, Particle* p2, Particle* p3, real l)
      : Super(p1, p2, p3, l) {}
  ParticleTripleLambda(Particle &p1, Particle &p2, Particle &p3, real &l)
      : Super(&p1, &p2, &p3, l) {}
};

class FixedTripleListLambda: public esutil::ESPPContainer< std::vector< ParticleTripleLambda > > {
 protected:
  typedef boost::unordered_multimap<longint, std::pair <longint, std::pair<longint, real> > > GlobalTriples;
  boost::signals2::connection sigAfterRecv, sigOnParticleChanged, sigBeforeSend;
  shared_ptr <storage::Storage> storage;
  GlobalTriples globalTriples;

 public:
  FixedTripleListLambda(shared_ptr <storage::Storage> _storage, real initLambda);
  ~FixedTripleListLambda();

  bool add(longint pid1, longint pid2, longint pid3);
  bool iadd(longint pid1, longint pid2, longint pid3);
  void beforeSendParticles(ParticleList &pl, class OutBuffer &buf);
  void afterRecvParticles(ParticleList &pl, class InBuffer &buf);
  void onParticlesChanged();
  bool remove(longint pid1, longint pid2, longint pid3);
  void updateParticlesStorage();

  python::list getTriples();
  python::list getTriplesLambda();

  real getLambda(longint pid1, longint pid2, longint pid3);
  void setLambda(longint pid1, longint pid2, longint pid3, real lambda);
  void setLambdaAll(real lambda);

    /** Get the number of triplets in the GlobalTriples list */
  int size() { return globalTriples.size(); }
  int totalSize();

  boost::signals2::signal3 <void, longint, longint, longint> onTupleAdded;
  boost::signals2::signal3 <void, longint, longint, longint> onTupleRemoved;

  static void registerPython();
 private:
  real initLambda_;
  static LOG4ESPP_DECL_LOGGER(theLogger);

};
}

#endif

