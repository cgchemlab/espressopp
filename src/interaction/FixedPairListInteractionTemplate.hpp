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
#ifndef _INTERACTION_FIXEDPAIRLISTINTERACTIONTEMPLATE_HPP
#define _INTERACTION_FIXEDPAIRLISTINTERACTIONTEMPLATE_HPP

#include "mpi.hpp"
#include "Interaction.hpp"
#include "Real3D.hpp"
#include "Tensor.hpp"
#include "Particle.hpp"
#include "FixedPairList.hpp"
#include "FixedPairListAdress.hpp"
#include "esutil/Array2D.hpp"
#include "bc/BC.hpp"
#include "SystemAccess.hpp"
#include "Interaction.hpp"
#include "types.hpp"

namespace espressopp {
  namespace interaction {
    template < typename _Potential >
    class FixedPairListInteractionTemplate: public Interaction, SystemAccess {
        
    protected:
      typedef _Potential Potential;
      
    public:
      FixedPairListInteractionTemplate
      (shared_ptr < System > system,
       shared_ptr < FixedPairList > _fixedpairList,
       shared_ptr < Potential > _potential)
        : SystemAccess(system), fixedpairList(_fixedpairList),
          potential(_potential)
      {
        if (! potential) {
          LOG4ESPP_ERROR(theLogger, "NULL potential");
        }
        precompute_pressure_ = false;
        precompute_energy_ = false;
      }

      virtual ~FixedPairListInteractionTemplate() {};

      void
      setFixedPairList(shared_ptr < FixedPairList > _fixedpairList) {
        fixedpairList = _fixedpairList;
      }

      shared_ptr < FixedPairList > getFixedPairList() {
        return fixedpairList;
      }

      void
      setPotential(shared_ptr < Potential> _potential) {
        if (_potential) {
          potential = _potential;
        } else {
          LOG4ESPP_ERROR(theLogger, "NULL potential");
        }
      }

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
      virtual int bondType() { return Pair; }

    protected:
      int ntypes;
      shared_ptr < FixedPairList > fixedpairList;
      shared_ptr < Potential > potential;
    };

    //////////////////////////////////////////////////
    // INLINE IMPLEMENTATION
    //////////////////////////////////////////////////
    template < typename _Potential > inline void
    FixedPairListInteractionTemplate < _Potential >::addForces() {
      LOG4ESPP_INFO(_Potential::theLogger, "adding forces of FixedPairList");
      const bc::BC& bc = *getSystemRef().bc;  // boundary conditions
      real ltMaxBondSqr = fixedpairList->getLongtimeMaxBondSqr();
      
      if (precompute_energy_) {
        e_local = 0.0;
      }
      if (precompute_pressure_) {
        w_virial = 0.0;
        w_wlocal = 0.0;
      }

      for (FixedPairList::PairList::Iterator it(*fixedpairList); it.isValid(); ++it) {
        Particle &p1 = *it->first;
        Particle &p2 = *it->second;
        Real3D dist;
        bc.getMinimumImageVector(dist, p1.position(), p2.position());
        Real3D force;
        real d = dist.sqr();
        if (d > ltMaxBondSqr) {
        	fixedpairList->setLongtimeMaxBondSqr(d);
        	ltMaxBondSqr = d;
        }
        if(potential->_computeForce(force, dist)) {
          /*if (force.isNaNInf()) {
            std::cout << " f: " << force
              << " p1 " << p1.id() << "p1.pos=" << p1.position()
              << " p2 " << p2.id() << " p2.pos=" << p2.position()
              << " d=" << d
              << std::endl;
            exit(1);
          }*/
          p1.force() += force;
          p2.force() -= force;
        }
        // Precompute local value of energy.
        if (precompute_energy_)
          e_local += potential->_computeEnergy(dist);
        // Precompute local value of 
        if (precompute_pressure_) {
          w_virial += dist * force;
          w_wlocal += Tensor(dist, force);
        }
      }
    }
    
    template < typename _Potential > inline real
    FixedPairListInteractionTemplate < _Potential >::computeEnergy() {

      LOG4ESPP_INFO(theLogger, "compute energy of the FixedPairList pairs");

      if (!precompute_energy_) {
        e_local = 0.0;
        const bc::BC& bc = *getSystemRef().bc;  // boundary conditions
        for (FixedPairList::PairList::Iterator it(*fixedpairList); it.isValid(); ++it) {
          const Particle &p1 = *it->first;
          const Particle &p2 = *it->second;
          Real3D r21;
          bc.getMinimumImageVectorBox(r21, p1.position(), p2.position());
          e_local += potential->_computeEnergy(r21);
        }
      }

      real esum;
      boost::mpi::all_reduce(*mpiWorld, e_local, esum, std::plus<real>());
      return esum;
    }
    
    template < typename _Potential > inline real
    FixedPairListInteractionTemplate < _Potential >::
    computeEnergyAA() {
      std::cout << "Warning! At the moment computeEnergyAA() in FixedPairListInteractionTemplate does not work." << std::endl;
      return 0.0;
    }
    
    template < typename _Potential > inline real
    FixedPairListInteractionTemplate < _Potential >::
    computeEnergyCG() {
      std::cout << "Warning! At the moment computeEnergyCG() in FixedPairListInteractionTemplate does not work." << std::endl;
      return 0.0;
    }   
    
    template < typename _Potential >
    inline void
    FixedPairListInteractionTemplate < _Potential >::
    computeVirialX(std::vector<real> &p_xx_total, int bins) {
              LOG4ESPP_INFO(theLogger, "compute virial p_xx of the pressure tensor slabwise");
    }

       
    template < typename _Potential > inline real
    FixedPairListInteractionTemplate < _Potential >::
    computeVirial() {
      LOG4ESPP_INFO(theLogger, "compute the virial for the FixedPair List");
     
      if (!precompute_pressure_) {
        w_virial = 0.0;
        const bc::BC& bc = *getSystemRef().bc;  // boundary conditions
        for (FixedPairList::PairList::Iterator it(*fixedpairList);
             it.isValid(); ++it) {                                         
          const Particle &p1 = *it->first;                                       
          const Particle &p2 = *it->second;                                      

          Real3D r21;
          bc.getMinimumImageVectorBox(r21, p1.position(), p2.position());
          Real3D force;
          if(potential->_computeForce(force, r21)) {
            w_virial += r21 * force;
          }
        }
      }
      
      real wsum;
      boost::mpi::all_reduce(*mpiWorld, w_virial, wsum, std::plus<real>());
      return wsum;
    }

    template < typename _Potential > inline void
    FixedPairListInteractionTemplate < _Potential >::computeVirialTensor(Tensor& w){
      LOG4ESPP_INFO(theLogger, "compute the virial tensor for the FixedPair List");
      
      if (!precompute_pressure_) {
        w_wlocal = 0.0;
        const bc::BC& bc = *getSystemRef().bc;  // boundary conditions
        for (FixedPairList::PairList::Iterator it(*fixedpairList);
             it.isValid(); ++it) {
          const Particle &p1 = *it->first;
          const Particle &p2 = *it->second;
          Real3D r21;
          bc.getMinimumImageVectorBox(r21, p1.position(), p2.position());
          Real3D force;
          if(potential->_computeForce(force, r21)) { 
            w_wlocal += Tensor(r21, force);
          }
        }
      }
      
      // reduce over all CPUs
      Tensor wsum(0.0);
      boost::mpi::all_reduce(*mpiWorld, (double*)&w_wlocal, 6, (double*)&wsum, std::plus<double>());
      w += wsum;
    }

    template < typename _Potential > inline void
    FixedPairListInteractionTemplate < _Potential >::
    computeVirialTensor(Tensor& w, real z){
      LOG4ESPP_INFO(theLogger, "compute the virial tensor for the FixedPair List");

      Tensor wlocal(0.0);
      const bc::BC& bc = *getSystemRef().bc;  // boundary conditions
      for (FixedPairList::PairList::Iterator it(*fixedpairList);
           it.isValid(); ++it) {
        const Particle &p1 = *it->first;
        const Particle &p2 = *it->second;
        Real3D p1pos = p1.position();
        Real3D p2pos = p2.position();
        
        if(  (p1pos[2]>=z && p2pos[2]<=z) ||
             (p1pos[2]<=z && p2pos[2]>=z) ){
          Real3D r21;
          bc.getMinimumImageVectorBox(r21, p1pos, p2pos);
          Real3D force;
          if(potential->_computeForce(force, r21)) { 
            wlocal += Tensor(r21, force);
          }
        }
      }
      
      // reduce over all CPUs
      Tensor wsum(0.0);
      boost::mpi::all_reduce(*mpiWorld, (double*)&wlocal, 6, (double*)&wsum, std::plus<double>());
      w += wsum;
    }
    
    template < typename _Potential > inline void
    FixedPairListInteractionTemplate < _Potential >::
    computeVirialTensor(Tensor *w, int n){
      LOG4ESPP_INFO(theLogger, "compute the virial tensor for the FixedPair List");

      const bc::BC& bc = *getSystemRef().bc;  // boundary conditions
      Real3D Li = bc.getBoxL();
      Tensor *wlocal = new Tensor[n];
      for(int i=0; i<n; i++) wlocal[i] = Tensor(0.0);
      for (FixedPairList::PairList::Iterator it(*fixedpairList);
           it.isValid(); ++it) {
        const Particle &p1 = *it->first;
        const Particle &p2 = *it->second;
        Real3D p1pos = p1.position();
        Real3D p2pos = p2.position();
        
        int position1 = (int)( n * p1pos[2]/Li[2]);
        int position2 = (int)( n * p1pos[2]/Li[2]);
        
        int maxpos = std::max(position1, position2);
        int minpos = std::min(position1, position2); 
        
        Real3D r21;
        bc.getMinimumImageVectorBox(r21, p1pos, p2pos);
        Real3D force;
        Tensor ww;
        if(potential->_computeForce(force, r21)) { 
          ww = Tensor(r21, force);
        }
        
        int i = minpos + 1;
        while(i<=maxpos){
          wlocal[i] += ww;
          i++;
        }
      }
      
      Tensor *wsum = new Tensor[n];
      boost::mpi::all_reduce(*mpiWorld, (double*)&wlocal, n, (double*)&wsum, std::plus<double>());
      
      for(int j=0; j<n; j++){
        w[j] += wsum[j];
      }

      delete [] wsum;
      delete [] wlocal;
    }
    
    template < typename _Potential >
    inline real
    FixedPairListInteractionTemplate< _Potential >::
    getMaxCutoff() {
      return potential->getCutoff();
    }
  }
}
#endif
