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

// ESPP_CLASS
#ifndef _INTERACTION_FIXEDTRIPLELISTINTERACTIONTEMPLATE_HPP
#define _INTERACTION_FIXEDTRIPLELISTINTERACTIONTEMPLATE_HPP

#include "mpi.hpp"
#include "Interaction.hpp"
#include "Real3D.hpp"
#include "Tensor.hpp"
#include "Particle.hpp"
#include "FixedTripleList.hpp"
#include "FixedTripleListAdress.hpp"
#include "esutil/Array2D.hpp"
#include "bc/BC.hpp"
#include "SystemAccess.hpp"
#include "types.hpp"

namespace espressopp {
  namespace interaction {
    template < typename _AngularPotential >
    class FixedTripleListInteractionTemplate : public Interaction, SystemAccess {
        
    protected:
      typedef _AngularPotential Potential;
      
    public:
      FixedTripleListInteractionTemplate
      (shared_ptr < System > _system,
       shared_ptr < FixedTripleList > _fixedtripleList,
       shared_ptr < Potential > _potential)
        : SystemAccess(_system), fixedtripleList(_fixedtripleList),
          potential(_potential)
      {
          if (! potential) {
                LOG4ESPP_ERROR(theLogger, "NULL potential");
          }
        //potentialArray = esutil::Array2D<Potential, esutil::enlarge>(0, 0, Potential());
      }

      virtual ~FixedTripleListInteractionTemplate() {};

      void
      setFixedTripleList(shared_ptr < FixedTripleList > _fixedtripleList) {
        fixedtripleList = _fixedtripleList;
      }

      shared_ptr < FixedTripleList > getFixedTripleList() {
        return fixedtripleList;
      }

      /*void
      setPotential(int type1, int type2, const Potential &potential) {
        potentialArray.at(type1, type2) = potential;
      }*/
      void
      setPotential(shared_ptr < Potential> _potential) {
         if (_potential) {
            potential = _potential;
         } else {
            LOG4ESPP_ERROR(theLogger, "NULL potential");
         }
      }

      /*Potential &getPotential(int type1, int type2) {
        return potentialArray.at(0, 0);
      }*/

      shared_ptr < Potential > getPotential() {
        return potential;
      }

      virtual void addForces();
      virtual real computeEnergy();
      virtual real computeEnergyAA();
      virtual real computeEnergyCG();      
      virtual void computeVirialX(std::vector<real> &p_xx_total, int bins); 
      virtual real computeVirial();
      virtual void computeVirialTensor(Tensor& w);
      virtual void computeVirialTensor(Tensor& w, real z);
      virtual void computeVirialTensor(Tensor *w, int n);
      virtual real getMaxCutoff();
      virtual int bondType() { return Angular; }

    protected:
      int ntypes;
      shared_ptr<FixedTripleList> fixedtripleList;
      //esutil::Array2D<Potential, esutil::enlarge> potentialArray;
      shared_ptr < Potential > potential;
    };

    //////////////////////////////////////////////////
    // INLINE IMPLEMENTATION
    //////////////////////////////////////////////////
    template < typename _AngularPotential > inline void
    FixedTripleListInteractionTemplate <_AngularPotential>::
    addForces() {
      LOG4ESPP_INFO(theLogger, "add forces computed by FixedTripleList");
      const bc::BC& bc = *getSystemRef().bc;  // boundary conditions
      if (precompute_energy_)
        e_local = 0.0;
      if (precompute_pressure_) {
        w_virial = 0.0;
        w_wlocal = 0.0;
      }
      for (FixedTripleList::TripleList::Iterator it(*fixedtripleList); it.isValid(); ++it) {
        Particle &p1 = *it->first;
        Particle &p2 = *it->second;
        Particle &p3 = *it->third;
        //const Potential &potential = getPotential(p1.type(), p2.type());
        Real3D dist12, dist32;
        bc.getMinimumImageVectorBox(dist12, p1.position(), p2.position());
        bc.getMinimumImageVectorBox(dist32, p3.position(), p2.position());
        Real3D force12, force32;
        potential->_computeForce(force12, force32, dist12, dist32);
        /*
        if (force12.isNaNInf() || force32.isNaNInf()) {
          std::cout << " f12: " << force12 << " f32: " << force32
          << " p1 " << p1.id() << "p1.pos=" << p1.position()
          << " p2 " << p2.id() << " p2.pos=" << p2.position()
          << " p3 " << p3.id() << " p2.pos=" << p3.position()
          << std::endl;
          exit(1);
        }*/

        p1.force() += force12;
        p2.force() -= force12 + force32;
        p3.force() += force32;
        if (precompute_energy_)
          e_local += potential->_computeEnergy(dist12, dist32);
        if (precompute_pressure_) {
          w_virial += dist12 * force12 + dist32 * force32;
          w_wlocal += Tensor(dist12, force12) + Tensor(dist32, force32);
        }
      }
    }

    template < typename _AngularPotential > inline real
    FixedTripleListInteractionTemplate < _AngularPotential >::
    computeEnergy() {
      LOG4ESPP_INFO(theLogger, "compute energy of the triples");
      
      if (!precompute_energy_) {
        const bc::BC& bc = *getSystemRef().bc;
        e_local = 0.0;
        for (FixedTripleList::TripleList::Iterator it(*fixedtripleList); it.isValid(); ++it) {
          const Particle &p1 = *it->first;
          const Particle &p2 = *it->second;
          const Particle &p3 = *it->third;
          //const Potential &potential = getPotential(p1.type(), p2.type());
          Real3D dist12 = bc.getMinimumImageVector(p1.position(), p2.position());
          Real3D dist32 = bc.getMinimumImageVector(p3.position(), p2.position());
          e_local += potential->_computeEnergy(dist12, dist32);
        }
      }
      real esum;
      boost::mpi::all_reduce(*mpiWorld, e_local, esum, std::plus<real>());
      return esum;
    }
    
    template < typename _AngularPotential > inline real
    FixedTripleListInteractionTemplate < _AngularPotential >::
    computeEnergyAA() {
      return 0.0;
    }
    
    template < typename _AngularPotential > inline real
    FixedTripleListInteractionTemplate < _AngularPotential >::
    computeEnergyCG() {
      return 0.0;
    }
           
    template < typename _AngularPotential >
    inline void
    FixedTripleListInteractionTemplate < _AngularPotential >::
    computeVirialX(std::vector<real> &p_xx_total, int bins) {
    }

    template < typename _AngularPotential > inline real
    FixedTripleListInteractionTemplate < _AngularPotential >::
    computeVirial() {
      LOG4ESPP_INFO(theLogger, "compute scalar virial of the triples");
      
      if (!precompute_pressure_) {
        const bc::BC& bc = *getSystemRef().bc;
        w_virial = 0.0;
        for (FixedTripleList::TripleList::Iterator it(*fixedtripleList); it.isValid(); ++it) {
          const Particle &p1 = *it->first;
          const Particle &p2 = *it->second;
          const Particle &p3 = *it->third;
          //const Potential &potential = getPotential(p1.type(), p2.type());
          const espressopp::bc::BC& bc = *getSystemRef().bc;
          Real3D dist12, dist32;
          bc.getMinimumImageVectorBox(dist12, p1.position(), p2.position());
          bc.getMinimumImageVectorBox(dist32, p3.position(), p2.position());
          Real3D force12, force32;
          potential->_computeForce(force12, force32, dist12, dist32);
          w_virial += dist12 * force12 + dist32 * force32;
        }
      }
      real wsum;
      boost::mpi::all_reduce(*mpiWorld, w_virial, wsum, std::plus<real>());
      return wsum;
    }

    template < typename _AngularPotential > inline void
    FixedTripleListInteractionTemplate < _AngularPotential >::
    computeVirialTensor(Tensor& w) {
      LOG4ESPP_INFO(theLogger, "compute the virial tensor of the triples");

      if (precompute_pressure_) {
        w_wlocal = 0.0;
        const bc::BC& bc = *getSystemRef().bc;
        for (FixedTripleList::TripleList::Iterator it(*fixedtripleList); it.isValid(); ++it){
          const Particle &p1 = *it->first;
          const Particle &p2 = *it->second;
          const Particle &p3 = *it->third;
          //const Potential &potential = getPotential(0, 0);
          Real3D r12, r32;
          bc.getMinimumImageVectorBox(r12, p1.position(), p2.position());
          bc.getMinimumImageVectorBox(r32, p3.position(), p2.position());
          Real3D force12, force32;
          potential->_computeForce(force12, force32, r12, r32);
          w_wlocal += Tensor(r12, force12) + Tensor(r32, force32);
        }
      }
      
      // reduce over all CPUs
      Tensor wsum(0.0);
      boost::mpi::all_reduce(*mpiWorld, (double*)&w_wlocal,6, (double*)&wsum, std::plus<double>());
      w += wsum;
    }

    template < typename _AngularPotential > inline void
    FixedTripleListInteractionTemplate < _AngularPotential >::
    computeVirialTensor(Tensor& w, real z) { }
    
    template < typename _AngularPotential > inline void
    FixedTripleListInteractionTemplate < _AngularPotential >::
    computeVirialTensor(Tensor *w, int n) { }
    
    template < typename _AngularPotential >
    inline real
    FixedTripleListInteractionTemplate< _AngularPotential >::
    getMaxCutoff() {
      return potential->getCutoff();
    }
  }
}
#endif
