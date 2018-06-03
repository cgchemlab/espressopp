/*
  Copyright (C) 2012,2013
      Max Planck Institute for Polymer Research
  Copyright (C) 2008,2009,2010,2011
      Max-Planck-Institute for Polymer Research & Fraunhofer SCAI
  
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
#include "ExtVelocity.hpp"

#include "types.hpp"
#include "System.hpp"
#include "storage/Storage.hpp"
#include "iterator/CellListIterator.hpp"

namespace espressopp {

  using namespace iterator;

  namespace integrator {

    LOG4ESPP_LOGGER(ExtVelocity::theLogger, "ExtVelocity");

    ExtVelocity::ExtVelocity(shared_ptr<System> system, const Real3D& _extVelocity)
    : Extension(system), extVelocity(_extVelocity)
    {
      LOG4ESPP_INFO(theLogger, "External Velocity for all particles constructed");
      allParticles = true;
    }

    ExtVelocity::ExtVelocity(shared_ptr<System> system, const Real3D& _extVelocity, shared_ptr< ParticleGroup > _particleGroup)
    : Extension(system), extVelocity(_extVelocity), particleGroup(_particleGroup)
    {
      LOG4ESPP_INFO(theLogger, "External Velocity for particle group constructed");
      allParticles = false;
    }

    void ExtVelocity::disconnect(){
      _aftInitF.disconnect();
    }

    void ExtVelocity::connect(){
      // connection to initialisation
      if (!allParticles) {
        _aftInitF  = integrator->befIntP.connect( boost::bind(&ExtVelocity::applyVelocityToGroup, this));
      } else {
    	_aftInitF  = integrator->befIntP.connect( boost::bind(&ExtVelocity::applyVelocityToAll, this));
      }
    }

    void ExtVelocity::setExtVelocity(Real3D& _extVelocity) {
    	extVelocity = _extVelocity;
    }

    Real3D& ExtVelocity::getExtVelocity() {
    	return extVelocity;
    }

    void ExtVelocity::setParticleGroup(shared_ptr< ParticleGroup > _particleGroup) {
        particleGroup = _particleGroup;
        if (allParticles) {
     	  disconnect();
     	  allParticles = false;
     	  connect();
        }
    }

     shared_ptr< ParticleGroup > ExtVelocity::getParticleGroup() {
    	//if (!allParticles) {
     	  return particleGroup;
    	//}
     }

     void ExtVelocity::applyVelocityToGroup() {
       LOG4ESPP_DEBUG(theLogger, "applying external force to particle group of size " << particleGroup->size());
       for (ParticleGroup::iterator it=particleGroup->begin(); it != particleGroup->end(); it++ ) {
    	 LOG4ESPP_DEBUG(theLogger, "applying external force to particle " << it->getId());
         it->velocity() += extVelocity;
       }
     }

     void ExtVelocity::applyVelocityToAll() {
       System& system = getSystemRef();
       CellList realCells = system.storage->getRealCells();
       for(CellListIterator cit(realCells); !cit.isDone(); ++cit) {
           cit->velocity() += extVelocity;
       }
     }

    /****************************************************
    ** REGISTRATION WITH PYTHON
    ****************************************************/

    void ExtVelocity::registerPython() {

      using namespace espressopp::python;

      class_<ExtVelocity, shared_ptr<ExtVelocity>, bases<Extension> >

        ("integrator_ExtVelocity", init< shared_ptr< System >, const Real3D& >())
        .def(init< shared_ptr< System >, const Real3D&, shared_ptr< ParticleGroup > >())
        .add_property("particleGroup", &ExtVelocity::getParticleGroup, &ExtVelocity::setParticleGroup)
        .def("getExtVelocity", &ExtVelocity::getExtVelocity, return_value_policy< reference_existing_object >())
        .def("setExtVelocity", &ExtVelocity::setExtVelocity )
        .def("connect", &ExtVelocity::connect)
        .def("disconnect", &ExtVelocity::disconnect)
        ;
    }

  }
}
