#ifndef _PAIRS_ALL_HPP
#define _PAIRS_ALL_HPP

#include "Set.hpp"
#include "Computer.hpp"
#include "Property.hpp"
#include "bc/BC.hpp"

namespace espresso {
  namespace pairs {
    /** Class that applies a Computer to all particle pairs of
	a given particle set.
    */
    class All : public Set {
    public:
      typedef shared_ptr< All > SelfPtr;

      static void registerPython();

      /** Constructor for this class 

	  \param bc are the boundary conditions that are needed for distance calculation.
	  \param set specifies the set of particles for which pairs will be considered.
	  \param coordinates the identifier of the coordinates property to use

      */
      All(bc::BC::SelfPtr _bc, 
	  particles::Set::SelfPtr _set, 
	  Property< Real3D >::SelfPtr _posProperty);

      using Set::foreach;

      /** This routine will apply a function operator to all pairs.
	  \param pairComputer is the object that provides the function to be applied to all pairs.
      */
      virtual void foreach(ApplyFunction _applyFunction);

      /** This routine will apply a function operator for read-only particles to all pairs.
	  \param pairComputer is the object that provides the read-only function to be applied to all pairs.
      */
      virtual void foreach(ConstApplyFunction _applyFunction) const;

    private:
      particles::Set::SelfPtr set;
      bc::BC::SelfPtr bc;
      Property< Real3D >::SelfPtr posProperty;
    };

  }
}

#endif
