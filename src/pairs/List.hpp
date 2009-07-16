#ifndef _PAIRS_LIST_HPP
#define _PAIRS_LIST_HPP

#include "Set.hpp"
#include "Computer.hpp"
#include "Particle.hpp"
#include "Property.hpp"
#include "bc/BC.hpp"

namespace espresso {
  namespace pairs {

    /** Class that applies a Computer to a list of pairs 
     */

    class List : public Set {
    public:
      typedef shared_ptr< List > SelfPtr;

    private:
      particles::Storage::SelfPtr storage; 
      bc::BC::SelfPtr bc;
      Property< Real3D >::SelfPtr coordinates;
      typedef std::pair< ParticleId,
			 ParticleId > Tuple;
      std::vector< Tuple > id_list;

    public:

      /** Destructor. */
      ~List();

      /** Constructor for this class 

	  \param bc are the boundary conditions that are needed for distance calculation.
	  \param storage specifies the particle storage to which particles belong
	  \param coordinates the identifier of the coordinates property to use

      */
      List (bc::BC::SelfPtr bc, 
            particles::Storage::SelfPtr storage,
	    Property< Real3D >::SelfPtr coordinates);

      size_t size() const;

      /** Ask if a particle pair tuple (id1, id2) is in the pair list

	  \param id1 is the identificaiton of particle 1
	  \param id2 is the identificaiton of particle 2
	  \return true if tuple (id1, id2) is in the list
 
      */

      bool findPair(ParticleId id1, ParticleId id2) const;

      /** Adding a particle pair tuple (id1, id2) to the pair list

	  \param id1 is the identificaiton of particle 1
	  \param id2 is the identificaiton of particle 2
 
	  Note: a tuple (id1, id2) can be added several times.
      */

      void addPair(ParticleId id1, ParticleId id2);

      /** Deleting a particle pair tuple (id1, id2) from the pair list

	  \param id1 is the identificaiton of particle 1
	  \param id2 is the identificaiton of particle 2
 
	  Particle (id1, i2) must be in the pair list otherwise exception.
      */

      void deletePair(ParticleId id1, ParticleId id2);

      /** This routine will apply a function operator to all pairs.

	  \param pairComputer is the object that provides the function to be applied to all pairs.

      */

      virtual void foreach(Computer& pairComputer);

      /** This routine will apply a function operator for read-only particles to all pairs.

	  \param pairComputer is the object that provides the read-only function to be applied to all pairs.

      */

      virtual void foreach(ConstComputer& pairComputer) const;

      /** Getter routine for the storage */

      particles::Storage::SelfPtr getStorage() const { return storage; }

    private:
      template< class Computer, class >
      void foreach(Computer& pairComputer) const;

    public:

      static void registerPython();

    };

  }
}

#endif
