/*
  Copyright (C) 2018
      Zidan Zhang (zidan.zhang@kuleuven.be)
  
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
#ifndef _INTEGRATOR_GeneralReactionScheme_HPP
#define _INTEGRATOR_GeneralReactionScheme_HPP

#include "Extension.hpp"
#include "boost/signals2.hpp"
#include "Real3D.hpp"
#include "Particle.hpp"
#include "iterator/CellListIterator.hpp"
#include "bc/BC.hpp"


namespace espressopp {
namespace integrator {

class GeneralReactionScheme: public Extension {
 public:
  /*** Select particles with probability p. */
  GeneralReactionScheme(shared_ptr<System> system, int interval);

  virtual ~GeneralReactionScheme() {
    disconnect();
  };

  /** Register this class so it can be used from Python. */
  static void registerPython();

 private:
  class ParticleData {
   public:
    ParticleData(Particle &p) {
      id = p.id();
      type = p.type();
      position = p.position();
      state = p.state();
      res_id = p.res_id();
    }
    ParticleData() {
      id = 0;
      type = 0;
      position = Real3D(0.0);
      state = 0;
      res_id = 0;
    }
    size_t id;
    longint type;
    Real3D position;
    longint state;
    longint res_id;
   private:
    friend class boost::serialization::access;
    template<class Archive>
    /**
     * Serialization to send data through MPI
     *
     * @tparam Archive
     * @param ar
     * @param version
     */
    void serialize(Archive & ar, const unsigned int version) {
      ar & id;
      ar & type;
      ar & position;
      ar & state;
      ar & res_id;
    }
  };

  boost::signals2::connection sig_aftIntV;
  int interval_;
  shared_ptr<bc::BC> bc_;

  void connect();
  void disconnect();

  /** Logger */
  static LOG4ESPP_DECL_LOGGER(theLogger);
  void updateParticles();
};

}
}

#endif
