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
#ifndef _FIXEDPAIRLAMBDALIST_HPP
#define _FIXEDPAIRLAMBDALIST_HPP

#include "log4espp.hpp"
#include "python.hpp"
#include "Particle.hpp"
#include "esutil/ESPPIterator.hpp"
#include <map>
#include <boost/signals2.hpp>
#include "types.hpp"

namespace espressopp {

class ParticlePairLambda : public Triple< class Particle*, class Particle*, real > {
 private:
  typedef Triple< class Particle*, class Particle*, real> Super;
 public:
  ParticlePairLambda() : Super() {}
  ParticlePairLambda(Particle* p1, Particle* p2, real l)
      : Super(p1, p2, l) {}
  ParticlePairLambda(Particle &p1, Particle &p2, real &l)
      : Super(&p1, &p2, l) {}
};

class FixedPairLambdaList: public esutil::ESPPContainer< std::vector< ParticlePairLambda > > {
 protected:
  typedef std::multimap <longint, std::pair<longint, real> > PairsLambda;
  boost::signals2::connection con1, con2, con3;
  shared_ptr <storage::Storage> storage;
  PairsLambda pairsLambda;
  real longtimeMaxBondSqr;

 public:
  FixedPairLambdaList(shared_ptr <storage::Storage> _storage, real initLambda);
  virtual ~FixedPairLambdaList();

  virtual bool add(longint pid1, longint pid2);
  virtual void beforeSendParticles(ParticleList &pl, class OutBuffer &buf);
  void afterRecvParticles(ParticleList &pl, class InBuffer &buf);
  virtual void onParticlesChanged();

  python::list getPairs();
  python::list getPairsLambda();

  real getLambda(longint pid1, longint pid2);
  /** Get the number of bonds in the GlobalPairs list */
  int size() { return pairsLambda.size(); }

  real getLongtimeMaxBondSqr() { return longtimeMaxBondSqr; }
  void setLongtimeMaxBondSqr(real d) { longtimeMaxBondSqr = d; };
  void resetLongtimeMaxBondSqr() { longtimeMaxBondSqr = 0.0; }

  boost::signals2::signal2 <void, longint, longint> onTupleAdded;

  static void registerPython();
 private:
  real initLambda_;
  static LOG4ESPP_DECL_LOGGER(theLogger);
};
}

#endif

