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
#ifndef _INTERACTION_FIXEDQUADRUPLELISTINTERACTIONTEMPLATE_HPP
#define _INTERACTION_FIXEDQUADRUPLELISTINTERACTIONTEMPLATE_HPP

#include "mpi.hpp"
#include "Interaction.hpp"
#include "Real3D.hpp"
#include "Tensor.hpp"
#include "Particle.hpp"
#include "FixedQuadrupleList.hpp"
#include "FixedQuadrupleListAdress.hpp"
#include "esutil/Array2D.hpp"
#include "bc/BC.hpp"
#include "SystemAccess.hpp"
#include "types.hpp"

namespace espressopp {
  namespace interaction {
    template < typename _DihedralPotential >
    class FixedQuadrupleListInteractionTemplate : public Interaction, SystemAccess {
    protected:
      typedef _DihedralPotential Potential;
    public:
      FixedQuadrupleListInteractionTemplate
      (shared_ptr < System > _system,
       shared_ptr < FixedQuadrupleList > _fixedquadrupleList,
       shared_ptr < Potential > _potential)
        : SystemAccess(_system), fixedquadrupleList(_fixedquadrupleList),
          potential(_potential)
      {
          if (! potential) {
              LOG4ESPP_ERROR(theLogger, "NULL potential");
          }

      }

      void
      setFixedQuadrupleList(shared_ptr < FixedQuadrupleList > _fixedquadrupleList) {
        fixedquadrupleList = _fixedquadrupleList;
      }

      virtual ~FixedQuadrupleListInteractionTemplate() {};

      shared_ptr < FixedQuadrupleList > getFixedQuadrupleList() {
        return fixedquadrupleList;
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
      virtual int bondType() { return Dihedral; }

    protected:
      int ntypes;
      shared_ptr < FixedQuadrupleList > fixedquadrupleList;
      shared_ptr < Potential > potential;
    };

    //////////////////////////////////////////////////
    // INLINE IMPLEMENTATION
    //////////////////////////////////////////////////
    template < typename _DihedralPotential > inline void
    FixedQuadrupleListInteractionTemplate < _DihedralPotential >::
    addForces() {

      LOG4ESPP_INFO(theLogger, "add forces computed by FixedQuadrupleList");

      const bc::BC& bc = *getSystemRef().bc;  // boundary conditions
      if (precompute_energy_)
        e_local = 0.0;
      if (precompute_pressure_) {
        w_virial = 0.0;
        w_wlocal = 0.0;
      }

      for (FixedQuadrupleList::QuadrupleList::Iterator it(*fixedquadrupleList); it.isValid(); ++it) {
        Particle &p1 = *it->first;
        Particle &p2 = *it->second;
        Particle &p3 = *it->third;
        Particle &p4 = *it->fourth;

        Real3D dist21, dist32, dist43; // 

        bc.getMinimumImageVectorBox(dist21, p2.position(), p1.position());
        bc.getMinimumImageVectorBox(dist32, p3.position(), p2.position());
        bc.getMinimumImageVectorBox(dist43, p4.position(), p3.position());

        Real3D force1, force2, force3, force4;  // result forces

        potential->_computeForce(force1, force2, force3, force4,
                                 dist21, dist32, dist43);
        /*
        if (force1.isNaNInf() || force2.isNaNInf() || force3.isNaNInf() || force4.isNaNInf()) {
          std::cout << " f1: " << force1 << " f2: " << force2
            << " f3:" << force3 << " f4:" << force4
            << " p1 " << p1.id() << "p1.pos=" << p1.position()
            << " p2 " << p2.id() << " p2.pos=" << p2.position()
            << " p3 " << p3.id() << " p3.pos=" << p3.position()
            << " p4 " << p4.id() << " p4.pos=" << p4.position()
            << std::endl;
          exit(1);
        }*/
        p1.force() += force1;
        p2.force() += force2; //p2.force() -= force2;
        p3.force() += force3;
        p4.force() += force4;

        if (precompute_energy_)
          e_local += potential->_computeEnergy(dist21, dist32, dist43);
        if (precompute_pressure_) {
          w_virial += dist21 * force1 + dist32 * force2;
          w_wlocal += Tensor(dist21, force1) - Tensor(dist32, force2);
        }
      }
    }

    template < typename _DihedralPotential >
    inline real
    FixedQuadrupleListInteractionTemplate < _DihedralPotential >::
    computeEnergy() {
      LOG4ESPP_INFO(theLogger, "compute energy of the quadruples");

      if (!precompute_energy_) {
        const bc::BC& bc = *getSystemRef().bc;  // boundary conditions
        e_local = 0.0;
        for (FixedQuadrupleList::QuadrupleList::Iterator it(*fixedquadrupleList); 
             it.isValid(); ++it) {
          const Particle &p1 = *it->first;
          const Particle &p2 = *it->second;
          const Particle &p3 = *it->third;
          const Particle &p4 = *it->fourth;

          Real3D dist21, dist32, dist43; // 

          bc.getMinimumImageVectorBox(dist21, p2.position(), p1.position());
          bc.getMinimumImageVectorBox(dist32, p3.position(), p2.position());
          bc.getMinimumImageVectorBox(dist43, p4.position(), p3.position());

          e_local += potential->_computeEnergy(dist21, dist32, dist43);
        }
      }
      real esum;
      boost::mpi::all_reduce(*mpiWorld, e_local, esum, std::plus<real>());
      return esum;
    }
    
    template < typename _DihedralPotential > inline real
    FixedQuadrupleListInteractionTemplate < _DihedralPotential >::
    computeEnergyAA() {
      return 0.0;
    }
    
    template < typename _DihedralPotential > inline real
    FixedQuadrupleListInteractionTemplate < _DihedralPotential >::
    computeEnergyCG() {
      return 0.0;
    }
    
    template < typename _DihedralPotential >
    inline void
    FixedQuadrupleListInteractionTemplate < _DihedralPotential >::
    computeVirialX(std::vector<real> &p_xx_total, int bins) {
    }

    template < typename _DihedralPotential >
    inline real
    FixedQuadrupleListInteractionTemplate < _DihedralPotential >::
    computeVirial() {
      LOG4ESPP_INFO(theLogger, "compute scalar virial of the quadruples");

      if (!precompute_pressure_) {
        w_virial = 0.0;
        const bc::BC& bc = *getSystemRef().bc;  // boundary conditions
        for (FixedQuadrupleList::QuadrupleList::Iterator it(*fixedquadrupleList); 
             it.isValid(); ++it) {
          const Particle &p1 = *it->first;
          const Particle &p2 = *it->second;
          const Particle &p3 = *it->third;
          const Particle &p4 = *it->fourth;

          Real3D dist21, dist32, dist43; 

          bc.getMinimumImageVectorBox(dist21, p2.position(), p1.position());
          bc.getMinimumImageVectorBox(dist32, p3.position(), p2.position());
          bc.getMinimumImageVectorBox(dist43, p4.position(), p3.position());

          Real3D force1, force2, force3, force4;

          potential->_computeForce(force1, force2, force3, force4,
                                  dist21, dist32, dist43);

          // TODO: formulas are not correct yet?

          w_virial += dist21 * force1 + dist32 * force2;
        }
      }
      
      real wsum = 0.0;
      boost::mpi::all_reduce(*mpiWorld, w_virial, wsum, std::plus<real>());
      return wsum;
    }

    template < typename _DihedralPotential >
    inline void
    FixedQuadrupleListInteractionTemplate < _DihedralPotential >::
    computeVirialTensor(Tensor& w) {
      LOG4ESPP_INFO(theLogger, "compute the virial tensor of the quadruples");
    
      if (precompute_pressure_) {
        w_wlocal = 0.0;
        const bc::BC& bc = *getSystemRef().bc;

        for (FixedQuadrupleList::QuadrupleList::Iterator it(*fixedquadrupleList); it.isValid(); ++it) {
          const Particle &p1 = *it->first;
          const Particle &p2 = *it->second;
          const Particle &p3 = *it->third;
          const Particle &p4 = *it->fourth;

          Real3D dist21, dist32, dist43; 

          bc.getMinimumImageVectorBox(dist21, p2.position(), p1.position());
          bc.getMinimumImageVectorBox(dist32, p3.position(), p2.position());
          bc.getMinimumImageVectorBox(dist43, p4.position(), p3.position());

          Real3D force1, force2, force3, force4;

          potential->_computeForce(force1, force2, force3, force4,
                                  dist21, dist32, dist43);

          // TODO: formulas are not correct yet

          w_wlocal += Tensor(dist21, force1) - Tensor(dist32, force2);
        }
      }
      // reduce over all CPUs
      Tensor wsum(0.0);
      boost::mpi::all_reduce(*mpiWorld, (double*)&w_wlocal, 6, (double*)&wsum, std::plus<double>());
      w += wsum;
    }


    template < typename _DihedralPotential >
    inline void
    FixedQuadrupleListInteractionTemplate < _DihedralPotential >::
    computeVirialTensor(Tensor& w, real z) {
      LOG4ESPP_INFO(theLogger, "compute the virial tensor of the quadruples");
    
      Tensor wlocal(0.0);
      const bc::BC& bc = *getSystemRef().bc;
      
      std::cout<<"Warning!!! computeVirialTensor in specified volume doesn't work for "
              "FixedQuadrupleListInteractionTemplate at the moment"<<std::endl;

      for (FixedQuadrupleList::QuadrupleList::Iterator it(*fixedquadrupleList); it.isValid(); ++it) {
        const Particle &p1 = *it->first;
        const Particle &p2 = *it->second;
        const Particle &p3 = *it->third;
        const Particle &p4 = *it->fourth;

        Real3D dist21, dist32, dist43; 

        bc.getMinimumImageVectorBox(dist21, p2.position(), p1.position());
        bc.getMinimumImageVectorBox(dist32, p3.position(), p2.position());
        bc.getMinimumImageVectorBox(dist43, p4.position(), p3.position());

        Real3D force1, force2, force3, force4;

        potential->_computeForce(force1, force2, force3, force4,
                                dist21, dist32, dist43);

        // TODO: formulas are not correct yet

        wlocal += Tensor(dist21, force1) - Tensor(dist32, force2);
      }
      // reduce over all CPUs
      Tensor wsum(0.0);
      boost::mpi::all_reduce(*mpiWorld, (double*)&wlocal, 6, (double*)&wsum, std::plus<double>());
      w += wsum;
    }

    template < typename _DihedralPotential >
    inline void
    FixedQuadrupleListInteractionTemplate < _DihedralPotential >::
    computeVirialTensor(Tensor *w, int n) {
      std::cout<<"Warning!!! computeVirialTensor in specified volume doesn't work for "
              "FixedQuadrupleListInteractionTemplate at the moment"<<std::endl;
    }
    
    template < typename _DihedralPotential >
    inline real
    FixedQuadrupleListInteractionTemplate< _DihedralPotential >::
    getMaxCutoff() {
      return potential->getCutoff();
    }
  }
}
#endif
