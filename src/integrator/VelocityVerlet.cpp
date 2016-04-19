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

//#include <iomanip>
#include "python.hpp"
#include "VelocityVerlet.hpp"
#include <iomanip>
#include "iterator/CellListIterator.hpp"
#include "interaction/Interaction.hpp"
#include "interaction/Potential.hpp"
#include "System.hpp"
#include "storage/Storage.hpp"
#include "mpi.hpp"

#ifdef VTRACE
#include "vampirtrace/vt_user.h"
#else
#define VT_TRACER( name)
#endif

namespace espressopp {
  using namespace std;
  namespace integrator {
    using namespace interaction;
    using namespace iterator;
    using namespace esutil;

    LOG4ESPP_LOGGER(VelocityVerlet::theLogger, "VelocityVerlet");

    VelocityVerlet::VelocityVerlet(shared_ptr< System > system) : MDIntegrator(system)
    {
      LOG4ESPP_INFO(theLogger, "construct VelocityVerlet");
      resortFlag = true;
      maxDist    = 0.0;
      timeIntegrate.reset();
      resetTimers();
      System& sys = getSystemRef();
    }

    VelocityVerlet::~VelocityVerlet()
    {
      LOG4ESPP_INFO(theLogger, "free VelocityVerlet");
    }

    void VelocityVerlet::run(int nsteps)
    {
      VT_TRACER("run");
      int nResorts = 0;
      real time;
      timeIntegrate.reset();
      System& system = getSystemRef();
      storage::Storage& storage = *system.storage;
      real skinHalf = 0.5 * system.getSkin();

      // Prepare the force comp timers if the size is not valid.
      const InteractionList& srIL = system.shortRangeInteractions;
      if (timeForceComp.size() < srIL.size()) {
        LOG4ESPP_DEBUG(theLogger, "Prepare timeForceComp");
        timeForceComp.clear();
        for (size_t i = 0; i < srIL.size(); i++) {
          timeForceComp.push_back(0.0);
        }
      }

      time = timeIntegrate.getElapsedTime();
      // signal
      runInit();
      timeRunInitS += timeIntegrate.getElapsedTime() - time;

      // Before start make sure that particles are on the right processor
      if (resortFlag) {
        VT_TRACER("resort");
        time = timeIntegrate.getElapsedTime();
        LOG4ESPP_INFO(theLogger, "resort particles");
        storage.decompose();
        maxDist = 0.0;
        resortFlag = false;
        timeResort += timeIntegrate.getElapsedTime();
      }

      bool recalcForces = true;  // TODO: more intelligent

      if (recalcForces) {
        LOG4ESPP_INFO(theLogger, "recalc forces before starting main integration loop");

        time = timeIntegrate.getElapsedTime();
        // signal
        recalc1();
        timeRecalc1S += timeIntegrate.getElapsedTime() - time;

        updateForces();
        if (LOG4ESPP_DEBUG_ON(theLogger)) {
            // printForces(false);   // forces are reduced to real particles
        }

        time = timeIntegrate.getElapsedTime();
        // signal
        recalc2();
        timeRecalc2S += timeIntegrate.getElapsedTime() - time;
      }

      LOG4ESPP_INFO(theLogger, "starting main integration loop (nsteps=" << nsteps << ")");
      
      for (int i = 0; i < nsteps; i++) {
        LOG4ESPP_INFO(theLogger, "Next step " << i << " of " << nsteps << " starts");

        //saveOldPos(); // save particle positions needed for constraints
        time = timeIntegrate.getElapsedTime();
        // signal
        befIntP();
        timeBefIntPS += timeIntegrate.getElapsedTime() - time;

        LOG4ESPP_INFO(theLogger, "updating positions and velocities")
        maxDist += integrate1();
        timeInt1 += timeIntegrate.getElapsedTime() - time;

        /*real cellsize = 1.4411685442;
        if (maxDist > 1.4411685442){
          cout<<"WARNING!!!!!! huge jump: "<<maxDist<<endl;
          exit(1);
        }*/

        time = timeIntegrate.getElapsedTime();
        // signal
        aftIntP();
        timeAftIntPS += timeIntegrate.getElapsedTime() - time;

        LOG4ESPP_INFO(theLogger, "maxDist = " << maxDist << ", skin/2 = " << skinHalf);

        if (maxDist > skinHalf) resortFlag = true;

        if (resortFlag) {
            VT_TRACER("resort1");
            time = timeIntegrate.getElapsedTime();
            LOG4ESPP_INFO(theLogger, "step " << i << ": resort particles");
            storage.decompose();
            maxDist  = 0.0;
            resortFlag = false;
            nResorts ++;
            timeResort += timeIntegrate.getElapsedTime() - time;
        }

        LOG4ESPP_INFO(theLogger, "updating forces")
        updateForces();

        timeIntegrate.startMeasure();
        // signal
        befIntV();
        timeBefIntVS += timeIntegrate.stopMeasure();

        time = timeIntegrate.getElapsedTime();
        integrate2();
        timeInt2 += timeIntegrate.getElapsedTime() - time;

        timeIntegrate.startMeasure();
        // signal
        aftIntV();
        aftIntV2();
        timeAftIntVS += timeIntegrate.stopMeasure();
      }

      timeRun = timeIntegrate.getElapsedTime();
      LOG4ESPP_INFO(theLogger, "finished run");
    }

    void VelocityVerlet::resetTimers() {
      // Prepare the force comp timers if the size is not valid.
      System& system = getSystemRef();
      const InteractionList& srIL = system.shortRangeInteractions;
      if (timeForceComp.size() < srIL.size()) {
        LOG4ESPP_DEBUG(theLogger, "Prepare timeForceComp");
        timeForceComp.clear();
        for (size_t i = 0; i < srIL.size(); i++) {
          timeForceComp.push_back(0.0);
        }
      }

      timeComm1  = 0.0;
      timeComm2  = 0.0;
      timeInt1   = 0.0;
      timeInt2   = 0.0;
      timeResort = 0.0;

      // Reset signal timers.
      timeRunInitS = 0.0;
      timeRecalc1S = 0.0;
      timeRecalc2S = 0.0;
      timeBefIntPS = 0.0;
      timeAftIntPS = 0.0;
      timeAftInitFS = 0.0;
      timeAftCalcFS = 0.0;
      timeBefIntVS = 0.0;
      timeAftIntVS = 0.0;
    }


    void VelocityVerlet::loadTimers(std::vector<real> &return_vector, std::vector<std::string> &labels) {

      return_vector.push_back(timeRun);
      labels.push_back("timeRun");
      for (int i = 0; i < timeForceComp.size(); i++) {
        std::stringstream ss;
        ss << "f" << i;
        return_vector.push_back(timeForceComp[i]);
        labels.push_back(ss.str());
      }

      // signal timers.
      return_vector.push_back(timeRunInitS);
      return_vector.push_back(timeRecalc1S);
      return_vector.push_back(timeRecalc2S);
      return_vector.push_back(timeBefIntPS);
      return_vector.push_back(timeAftIntPS);
      return_vector.push_back(timeAftCalcFS);
      return_vector.push_back(timeBefIntVS);
      return_vector.push_back(timeAftIntVS);

      return_vector.push_back(timeComm1);
      return_vector.push_back(timeComm2);
      return_vector.push_back(timeInt1);
      return_vector.push_back(timeInt2);
      return_vector.push_back(timeResort);

      labels.push_back("timeRunInitS");
      labels.push_back("timeRecalc1S");
      labels.push_back("timeRecalc2S");
      labels.push_back("timeBefIntPS");
      labels.push_back("timeAftIntPS");
      labels.push_back("timeAftCalcFS");
      labels.push_back("timeBefIntVS");
      labels.push_back("timeAftIntVS");

      labels.push_back("timeComm1");
      labels.push_back("timeComm2");
      labels.push_back("timeInt1");
      labels.push_back("timeInt2");
      labels.push_back("timeResort");
    }

    static boost::python::object wrapGetTimers(class VelocityVerlet* obj) {
      std::vector<real> timers;
      std::vector<std::string> labels;
      obj->loadTimers(timers, labels);

      boost::python::list return_list;
      for (int i = 0; i < timers.size(); i++) {
        return_list.append(boost::python::make_tuple(labels[i], timers[i]));
      }
      return return_list;
    }


    real VelocityVerlet::integrate1()
    {
      System& system = getSystemRef();
      CellList realCells = system.storage->getRealCells();

      // loop over all particles of the local cells
      int count = 0;
      real maxSqDist = 0.0; // maximal square distance a particle moves
      for(CellListIterator cit(realCells); !cit.isDone(); ++cit) {
        real sqDist = 0.0;
        LOG4ESPP_INFO(theLogger, "updating first half step of velocities and full step of positions")
        LOG4ESPP_DEBUG(theLogger, "Particle " << cit->id() <<
                ", pos = " << cit->position() <<
                ", v = " << cit->velocity() <<
                ", f = " << cit->force());

        /* more precise for DEBUG:
        printf("Particle %d, pos = %16.12f %16.12f %16.12f, v = %16.12f, %16.12f %16.12f, f = %16.12f %16.12f %16.12f\n",
            cit->p.id, cit->r.p[0], cit->r.p[1], cit->r.p[2],
                cit->m.v[0], cit->m.v[1], cit->m.v[2],
            cit->f.f[0], cit->f.f[1], cit->f.f[2]);
        */

        real dtfm = 0.5 * dt / cit->mass();

        // Propagate velocities: v(t+0.5*dt) = v(t) + 0.5*dt * f(t)
        cit->velocity() += dtfm * cit->force();

        // Propagate positions (only NVT): p(t + dt) = p(t) + dt * v(t+0.5*dt) 
        Real3D deltaP = cit->velocity();
        
        deltaP *= dt;
        cit->position() += deltaP;
        sqDist += deltaP * deltaP;

        count++;

        maxSqDist = std::max(maxSqDist, sqDist);
      }
      
      // signal
      inIntP(maxSqDist);

      real maxAllSqDist;
      mpi::all_reduce(*system.comm, maxSqDist, maxAllSqDist, boost::mpi::maximum<real>());

      LOG4ESPP_INFO(theLogger, "moved " << count << " particles in integrate1" <<
		    ", max move local = " << sqrt(maxSqDist) <<
		    ", global = " << sqrt(maxAllSqDist));
      
      return sqrt(maxAllSqDist);
    }

    void VelocityVerlet::integrate2()
    {
      LOG4ESPP_INFO(theLogger, "updating second half step of velocities")
      System& system = getSystemRef();
      CellList realCells = system.storage->getRealCells();

      // loop over all particles of the local cells
      real half_dt = 0.5 * dt; 
      for(CellListIterator cit(realCells); !cit.isDone(); ++cit) {
        real dtfm = half_dt / cit->mass();
        /* Propagate velocities: v(t+0.5*dt) = v(t) + 0.5*dt * f(t) */
        cit->velocity() += dtfm * cit->force();
      }
      
      step++;
    }

    void VelocityVerlet::calcForces()
    {
      VT_TRACER("forces");

      LOG4ESPP_INFO(theLogger, "calculate forces");

      initForces();

      timeIntegrate.startMeasure();
      // signal
      aftInitF();
      timeAftInitFS += timeIntegrate.stopMeasure();

      System& sys = getSystemRef();
      const InteractionList& srIL = sys.shortRangeInteractions;
      real time;
      for (size_t i = 0; i < srIL.size(); i++) {
        LOG4ESPP_INFO(theLogger, "compute forces for srIL " << i << " of " << srIL.size());
        time = timeIntegrate.getElapsedTime();
        srIL[i]->addForces();
        timeForceComp[i] += timeIntegrate.getElapsedTime() - time;
      }
    }

    void VelocityVerlet::updateForces()
    {
      LOG4ESPP_INFO(theLogger, "update ghosts, calculate forces and collect ghost forces")
      real time;
      storage::Storage& storage = *getSystemRef().storage;
      time = timeIntegrate.getElapsedTime();
      {
        VT_TRACER("commF");
        storage.updateGhosts();
      }
      timeComm1 += timeIntegrate.getElapsedTime() - time;
      time = timeIntegrate.getElapsedTime();
      calcForces();
      timeForce += timeIntegrate.getElapsedTime() - time;
      time = timeIntegrate.getElapsedTime();
      {
        VT_TRACER("commR");
        storage.collectGhostForces();
      }
      timeComm2 += timeIntegrate.getElapsedTime() - time;

      timeIntegrate.startMeasure();
      // signal
      aftCalcF();
      timeAftCalcFS += timeIntegrate.stopMeasure();
    }

    void VelocityVerlet::initForces()
    {
      // forces are initialized for real + ghost particles

      System& system = getSystemRef();
      CellList localCells = system.storage->getLocalCells();

      LOG4ESPP_INFO(theLogger, "init forces for real + ghost particles");

      for(CellListIterator cit(localCells); !cit.isDone(); ++cit) {
        cit->force() = 0.0;
        cit->drift() = 0.0;   // Can in principle be commented, when drift is not used.
      }
    }

    void VelocityVerlet::printForces(bool withGhosts)
    {
      // print forces of real + ghost particles

      System& system = getSystemRef();
      CellList cells;

      if (withGhosts) {
	    cells = system.storage->getLocalCells();
	    LOG4ESPP_DEBUG(theLogger, "local forces");
      } else {
	    cells = system.storage->getRealCells();
	    LOG4ESPP_DEBUG(theLogger, "real forces");
      }
  
      for(CellListIterator cit(cells); !cit.isDone(); ++cit) {
	    LOG4ESPP_DEBUG(theLogger, "Particle " << cit->id() << ", force = " << cit->force());
      }
    }

    void VelocityVerlet::printPositions(bool withGhosts)
    {
      // print positions of real + ghost particles

      System& system = getSystemRef();
      CellList cells;

      if (withGhosts) {
	    cells = system.storage->getLocalCells();
	    LOG4ESPP_DEBUG(theLogger, "local positions");
      } else {
	    cells = system.storage->getRealCells();
	    LOG4ESPP_DEBUG(theLogger, "real positions");
      }
  
      for(CellListIterator cit(cells); !cit.isDone(); ++cit) {
	    LOG4ESPP_DEBUG(theLogger, "Particle " << cit->id() << ", position = " << cit->position());
      }
    }

    /****************************************************
    ** REGISTRATION WITH PYTHON
    ****************************************************/

    void VelocityVerlet::registerPython() {

      using namespace espressopp::python;

      // Note: use noncopyable and no_init for abstract classes
      class_<VelocityVerlet, bases<MDIntegrator>, boost::noncopyable >
        ("integrator_VelocityVerlet", init< shared_ptr<System> >())
        .def("getTimers", &wrapGetTimers)
        .def("resetTimers", &VelocityVerlet::resetTimers)
        ;
    }
  }
}
