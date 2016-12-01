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

// ESPP_CLASS
#ifndef _FIXEDQUADRUPLELISTLAMBDA_HPP
#define _FIXEDQUADRUPLELISTLAMBDA_HPP

#include "log4espp.hpp"
#include "types.hpp"
#include "Triple.hpp"

#include "Particle.hpp"
#include "FixedQuadrupleList.hpp"
#include "esutil/ESPPIterator.hpp"
#include <boost/unordered_map.hpp>
#include <boost/signals2.hpp>
#include <boost/python/list.hpp>

namespace espressopp {

class ParticleQuadrupleLambda {
 public:
  ParticleQuadrupleLambda(Particle *p1, Particle *p2, Particle *p3, Particle *p4, real l)
      : first(p1), second(p2), third(p3), fourth(p4), lambda(l) {}

  ParticleQuadrupleLambda(Particle &p1, Particle &p2, Particle &p3, Particle &p4, real l)
      : first(&p1), second(&p2), third(&p3), fourth(&p4), lambda(l) {}

  Particle *first;
  Particle *second;
  Particle *third;
  Particle *fourth;
  real lambda;
};

class FixedQuadrupleListLambda : public FixedQuadrupleList {
 protected:
  typedef boost::unordered_multimap<longint, std::pair<Triple<longint, longint, longint>, real> > QuadruplesLambda;
  shared_ptr <storage::Storage> storage;
  using QuadrupleList::add;

 public:
  typedef std::vector<ParticleQuadrupleLambda> ParticleQuadruplesLambda;
  typedef esutil::ESPPIterator<ParticleQuadruplesLambda> IteratorParticleLambda;

 public:
  FixedQuadrupleListLambda(shared_ptr <storage::Storage> _storage, real lambda0);

  bool add(longint pid1, longint pid2, longint pid3, longint pid4);
  /// Non-blocking add method.
  bool iadd(longint pid1, longint pid2, longint pid3, longint pid4);

  bool remove(longint pid1, longint pid2, longint pid3, longint pid4);

  void beforeSendParticles(ParticleList &pl, class OutBuffer &buf);
  void afterRecvParticles(ParticleList &pl, class InBuffer &buf);
  virtual void onParticlesChanged();
  virtual void updateParticlesStorage();

  python::list getQuadruples();
  real getLambda(longint pid1, longint pid2, longint pid3, longint pid4);
  void setLambda(longint pid1, longint pid2, longint pid3, longint pid4, real lambda);
  void setAllLambda(real lambda);
  void incrementAllLambda(real d_lambda);

  ParticleQuadruplesLambda& getParticleQuadruples() { return particleQuadruplesLambda_; }

  boost::python::list getQuadruplesLambda();

  static void registerPython();

 private:
  real lambda0_;
  QuadruplesLambda quadruplesLambda_;
  ParticleQuadruplesLambda particleQuadruplesLambda_;
  static LOG4ESPP_DECL_LOGGER(theLogger);
  python::list getAllQuadruples();
};
}

#endif
