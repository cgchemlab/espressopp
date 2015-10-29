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
#ifndef _VERLETLIST_HPP
#define _VERLETLIST_HPP

#include "log4espp.hpp"
#include "types.hpp"
#include "python.hpp"
#include "Particle.hpp"
#include "Buffer.hpp"
#include "SystemAccess.hpp"
#include "integrator/MDIntegrator.hpp"
#include "boost/signals2.hpp"
#include "boost/unordered_set.hpp"

namespace espressopp {

class DynamicExcludeList {
public:
    typedef boost::unordered_set<std::pair<longint, longint> > ExcludeList;
    DynamicExcludeList(shared_ptr<integrator::MDIntegrator> integrator);
    ~DynamicExcludeList();
    void exclude(longint pid1, longint pid2);
    void unexclude(longint pid1, longint pid2);
    void connect();
    void disconnect();
    shared_ptr<ExcludeList> getExList { return exList; }

    static void registerPython();
private:
    shared_ptr<integrator::MDIntegrator> integrator_;
    shared_ptr<ExcludeList> exList;
    // Helper lists.
    ExcludeList exList_add;
    ExcludeList exList_remove;
    bool exListDirty;
    void updateList();

    boost::signals2::connection aftIntV;
    static LOG4ESPP_DECL_LOGGER(theLogger);
};

/** Class that builds and stores verlet lists.

    ToDo: register at system for rebuild

*/

  class VerletList : public SystemAccess {

  public:
    typedef boost::unordered_set<std::pair<longint, longint> > ExcludeList;
    /** Build a verlet list of all particle pairs in the storage
	whose distance is less than a given cutoff.

	\param system is the system for which the verlet list is built
	\param cut is the cutoff value for the 

    */

    VerletList(shared_ptr< System >, real cut, bool rebuildVL);
    VerletList(shared_ptr< System >, real cut, shared_ptr<DynamicExcludeList>, bool rebuildVL);

    ~VerletList();

    PairList& getPairs() { return vlPairs; }

    python::tuple getPair(int i);

    real getVerletCutoff(); // returns cutoff + skin

    void connect();

    void disconnect();

    void rebuild();

    /** Get the total number of pairs for the Verlet list */
    int totalSize() const;

    //** Get the number of pairs for the local Verlet list */
    int localSize() const;

    /** Add pairs to exclusion list */
    bool exclude(longint pid1, longint pid2);

    void unexclude(longint pid1, longint pid2);

    /** Get the number of times the Verlet list has been rebuilt */
    int getBuilds() const { return builds; }

    /** Set the number of times the Verlet list has been rebuilt */
    void setBuilds(int _builds) { builds = _builds; }

    /** Register this class so it can be used from Python. */
    static void registerPython();

  protected:

    void checkPair(Particle &pt1, Particle &pt2);
    void afterRecvParticles(ParticleList &unused_pl, InBuffer &unused_buf);
    PairList vlPairs;
    shared_ptr<ExcludeList> exList; // exclusion list
    bool dynamicExList;

    real cutsq;
    real cut;
    real cutVerlet;
    
    int builds;
    boost::signals2::connection connectionResort;

    static LOG4ESPP_DECL_LOGGER(theLogger);

  };

}

#endif
