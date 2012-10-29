// ESPP_CLASS
#ifndef _VERLETLISTTRIPLE_HPP
#define _VERLETLISTTRIPLE_HPP

#include "log4espp.hpp"
#include "types.hpp"
#include "python.hpp"
#include "Particle.hpp"
#include "SystemAccess.hpp"
#include "boost/signals2.hpp"
#include "boost/unordered_set.hpp"

namespace espresso {

/** Class that builds and stores verlet lists for 3-body interactions.
    ToDo: register at system for rebuild
*/

  class VerletListTriple : public SystemAccess {

  public:

    /** Build a verlet list of all particle triples in the storage
	whose distance is less than a given cutoff.

	\param system is the system for which the verlet list triples are built
	\param cut is the cutoff value for the 

    */

    VerletListTriple(shared_ptr< System >, real cut, bool rebuildVL);

    ~VerletListTriple();

    TripleList& getTriples() { return vlTriples; }

    python::tuple getTriple(int i);

    real getVerletCutoff(); // returns cutoff + skin

    void connect();

    void disconnect();

    void rebuild();

    /** Get the total number of triples for the Verlet Triple list */
    int totalSize() const;

    //** Get the number of triples for the local Verlet list */
    int localSize() const;
    
    /** Add Triples to exclusion list */
    bool exclude(longint pid1, longint pid2, longint pid3);

    /** Get the number of times the Verlet Triple list has been rebuilt */
    int getBuilds() const { return builds; }

    /** Set the number of times the Verlet Triple  list has been rebuilt */
    void setBuilds(int _builds) { builds = _builds; }

    /** Register this class so it can be used from Python. */
    static void registerPython();

  protected:

    void checkTriple(Particle &pt1, Particle &pt2, Particle &pt3);
    TripleList vlTriples;
    
    boost::unordered_set< longint> exList; // exclusion list
    
    real cutsq;
    real cut;
    real cutVerlet;
    
    int builds;
    boost::signals2::connection connectionResort;

    static LOG4ESPP_DECL_LOGGER(theLogger);
  };

}

#endif
