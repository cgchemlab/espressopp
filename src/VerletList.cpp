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
#include "VerletList.hpp"
#include "Real3D.hpp"
#include "Particle.hpp"
#include "Cell.hpp"
#include "System.hpp"
#include "storage/Storage.hpp"
#include "bc/BC.hpp"
#include "iterator/CellListAllPairsIterator.hpp"

namespace espressopp {
using namespace espressopp::iterator;

LOG4ESPP_LOGGER(DynamicExcludeList::theLogger, "VerletList");

DynamicExcludeList::DynamicExcludeList(shared_ptr<integrator::MDIntegrator> integrator):
    integrator_(integrator) {
  LOG4ESPP_INFO(theLogger, "construct of DynamicExcludeList");
  exListDirty = true;
  exList = boost::make_shared<ExcludeList>();
  connect();
}

DynamicExcludeList::~DynamicExcludeList() {
  disconnect();
}

void DynamicExcludeList::connect() {
  LOG4ESPP_INFO(theLogger, "Connected to integrator");
  aftIntV = integrator_->aftIntV.connect(boost::bind(&DynamicExcludeList::updateList, this));
}

void DynamicExcludeList::disconnect() {
  LOG4ESPP_INFO(theLogger, "Disconnected from integrator");
  aftIntV.disconnect();
}

void DynamicExcludeList::observe(shared_ptr<FixedPairList> fpl) {
  fpl->onTupleAdded.connect(
      boost::bind(&DynamicExcludeList::exclude, this, _1, _2));

  // TODO: Handle when bond is removed from observerd  fixed pair list. Currently this is not
  // implemented there.
}

void DynamicExcludeList::updateList() {
  LOG4ESPP_INFO(theLogger, "Update dynamic list.");
  // Collect state from all CPUs. If somewhere list is dirty then gather and scatter.
  bool global_exListDirty;
  mpi::all_reduce(*(integrator_->getSystem()->comm), exListDirty,
                  global_exListDirty, std::logical_or<bool>());
  if (global_exListDirty) {
    LOG4ESPP_INFO(theLogger, "Exclude dynamic list is dirty, exchange it.");
    std::vector<longint> out_buffer;
    std::vector<std::vector<longint> > in_buffer;

    // Prepare output.
    out_buffer.push_back(2*exList_remove.size());
    out_buffer.insert(out_buffer.end(), exList_remove.begin(), exList_remove.end());
    out_buffer.push_back(2*exList_add.size());
    out_buffer.insert(out_buffer.end(), exList_add.begin(), exList_add.end());
    // Gather everywhere everything.
    mpi::all_gather(*(integrator_->getSystem()->comm), out_buffer, in_buffer);
    //Update list.
    for (std::vector<std::vector<longint> >::iterator it = in_buffer.begin();
         it != in_buffer.end(); ++it) {
      for (int i = 1; i < (it->at(0)+1); i=i+2) {
        exList->erase(std::make_pair(it->at(i), it->at(i+1)));
        exList->erase(std::make_pair(it->at(i+1), it->at(i)));
        LOG4ESPP_DEBUG(theLogger, "removed pair: " << it->at(i) << "-" << it->at(i+1));
      }
      for (int i = (it->at(0)+2); i < it->size(); i=i+2) {
        LOG4ESPP_DEBUG(theLogger, "add pair: " << it->at(i) << "-" << it->at(i+1));
        exList->insert(std::make_pair(it->at(i), it->at(i+1)));
        exList->insert(std::make_pair(it->at(i+1), it->at(i)));
      }
    }
  }
  exList_remove.clear();
  exList_add.clear();
  exListDirty = false;
}
void DynamicExcludeList::setExListDirty(bool val) {
  if (!val) {
    exList_remove.clear();
    exList_add.clear();
  }
  exListDirty = val;
}

python::list DynamicExcludeList::getList() {
  python::list return_list;
  for (ExcludeList::iterator it = exList->begin(); it != exList->end(); ++it) {
    return_list.append(python::make_tuple(it->first, it->second));
  }
  return return_list;
}

void DynamicExcludeList::exclude(longint pid1, longint pid2) {
  exList->insert(std::make_pair(pid1, pid2));

  exList_add.push_back(pid1);
  exList_add.push_back(pid2);
  exListDirty = true;
}

void DynamicExcludeList::unexclude(longint pid1, longint pid2) {
  exList_remove.push_back(pid1);
  exList_remove.push_back(pid2);

  exList->erase(std::make_pair(pid1, pid2));
  exList->erase(std::make_pair(pid2, pid1));
  exListDirty = true;
}

void DynamicExcludeList::registerPython() {
  using namespace espressopp::python;

  class_<DynamicExcludeList, shared_ptr<DynamicExcludeList> >
      ("DynamicExcludeList", init< shared_ptr<integrator::MDIntegrator> >())
       .add_property("is_dirty", &DynamicExcludeList::getExListDirty,
                     &DynamicExcludeList::setExListDirty)
       .add_property("size", &DynamicExcludeList::getSize)
       .def("exclude", &DynamicExcludeList::exclude)
       .def("observe", &DynamicExcludeList::observe)
       .def("unexclude", &DynamicExcludeList::unexclude)
       .def("get_list", &DynamicExcludeList::getList)
       .def("connect", &DynamicExcludeList::connect)
       .def("disconnect", &DynamicExcludeList::disconnect);
}

/** Implementation of VerletList **/

  LOG4ESPP_LOGGER(VerletList::theLogger, "VerletList");
  // cut is a cutoff (without skin)
  VerletList::VerletList(shared_ptr<System> system, real _cut, bool rebuildVL) : SystemAccess(system)
  {
    LOG4ESPP_INFO(theLogger, "construct VerletList, cut = " << _cut);
  
    if (!system->storage) {
       throw std::runtime_error("system has no storage");
    }

    cut = _cut;
    cutVerlet = cut + system -> getSkin();
    cutsq = cutVerlet * cutVerlet;
    builds = 0;

    exList = boost::make_shared<ExcludeList>();
    isDynamicExList = false;

    if (rebuildVL) rebuild(); // not called if exclutions are provided

    // make a connection to System to invoke rebuild on resort
    connect();
  }
  
  VerletList::VerletList(shared_ptr<System> system, real _cut,
                         shared_ptr<DynamicExcludeList> dynamicExList_, bool rebuildVL):
      SystemAccess(system), dynamicExcludeList(dynamicExList_) {
    LOG4ESPP_INFO(theLogger, "construct VerletList with dynamic exclusion list, cut = " << _cut);

    if (!system->storage) {
      throw std::runtime_error("system has no storage");
    }

    cut = _cut;
    cutVerlet = cut + system -> getSkin();
    cutsq = cutVerlet * cutVerlet;
    builds = 0;

    exList = dynamicExList_->getExList();

    isDynamicExList = true;

    if (rebuildVL) rebuild(); // not called if exclutions are provided

    // make a connection to System to invoke rebuild on resort
    connect();
  }

  real VerletList::getVerletCutoff(){
    return cutVerlet;
  }
  
  void VerletList::connect()
  {
  // make a connection to System to invoke rebuild on resort
  connectionResort = getSystem()->storage->onParticlesChanged.connect(
      boost::bind(&VerletList::rebuild, this));
  }

  void VerletList::disconnect()
  {

  // disconnect from System to avoid rebuild on resort
  connectionResort.disconnect();
  }

  /*-------------------------------------------------------------*/
  
  void VerletList::rebuild()
  {
    //real cutVerlet = cut + getSystem() -> getSkin();
    cutVerlet = cut + getSystem() -> getSkin();
    cutsq = cutVerlet * cutVerlet;
    
    vlPairs.clear();

    // add particles to adress zone
    CellList cl = getSystem()->storage->getRealCells();
    LOG4ESPP_DEBUG(theLogger, "local cell list size = " << cl.size());
    for (CellListAllPairsIterator it(cl); it.isValid(); ++it) {
      checkPair(*it->first, *it->second);
      LOG4ESPP_DEBUG(theLogger, "checking particles " << it->first->id() << " and " << it->second->id());
    }
    
    builds++;
    LOG4ESPP_DEBUG(theLogger, "rebuilt VerletList (count=" << builds << "), cutsq = " << cutsq
                 << " local size = " << vlPairs.size());
  }
  

  /*-------------------------------------------------------------*/
  
  void VerletList::checkPair(Particle& pt1, Particle& pt2)
  {

    Real3D d = pt1.position() - pt2.position();
    real distsq = d.sqr();

    LOG4ESPP_TRACE(theLogger, "p1: " << pt1.id()
                   << " @ " << pt1.position() 
		   << " - p2: " << pt2.id() << " @ " << pt2.position()
		   << " -> distsq = " << distsq);

    if (distsq > cutsq) return;

    // see if it's in the exclusion list (both directions)
    if (exList->count(std::make_pair(pt1.id(), pt2.id())) == 1) return;
    if (exList->count(std::make_pair(pt2.id(), pt1.id())) == 1) return;

    vlPairs.add(pt1, pt2); // add pair to Verlet List
  }
  
  /*-------------------------------------------------------------*/
  
  int VerletList::totalSize() const
  {
    System& system = getSystemRef();
    int size = localSize();
    int allsize;
  
    mpi::all_reduce(*system.comm, size, allsize, std::plus<int>());
    return allsize;
  }

  int VerletList::localSize() const
  {
    System& system = getSystemRef();
    return vlPairs.size();
  }

  python::tuple VerletList::getPair(int i) {
	  if (i <= 0 || i > vlPairs.size()) {
	    std::cout << "ERROR VerletList pair " << i << " does not exists" << std::endl;
	    return python::make_tuple();
	  } else {
	    return python::make_tuple(vlPairs[i-1].first->id(), vlPairs[i-1].second->id());
	  }
  }


  bool VerletList::exclude(longint pid1, longint pid2) {
      if (isDynamicExList) {
        dynamicExcludeList->exclude(pid1, pid2);
      } else {
        exList->insert(std::make_pair(pid1, pid2));
      }
      return true;
  }

  bool VerletList::unexclude(longint pid1, longint pid2) {
    if (isDynamicExList) {
      dynamicExcludeList->unexclude(pid1, pid2);
    } else {
      exList->erase(std::make_pair(pid1, pid2));
      exList->erase(std::make_pair(pid2, pid1));
    }
  }
  

  /*-------------------------------------------------------------*/
  
  VerletList::~VerletList()
  {
    LOG4ESPP_INFO(theLogger, "~VerletList");
  
    if (!connectionResort.connected()) {
      connectionResort.disconnect();
    }
  }
  
  /****************************************************
  ** REGISTRATION WITH PYTHON
  ****************************************************/
  
  void VerletList::registerPython() {
    using namespace espressopp::python;

    bool (VerletList::*pyExclude)(longint pid1, longint pid2)
          = &VerletList::exclude;


    class_<VerletList, shared_ptr<VerletList> >
      ("VerletList", init< shared_ptr<System>, real, bool >())
      .def(init<shared_ptr<System>, real, shared_ptr<DynamicExcludeList>, bool>())
      .add_property("system", &SystemAccess::getSystem)
      .add_property("builds", &VerletList::getBuilds, &VerletList::setBuilds)
      .def("totalSize", &VerletList::totalSize)
      .def("localSize", &VerletList::localSize)
      .def("getPair", &VerletList::getPair)
      .def("exclude", pyExclude)
      .def("rebuild", &VerletList::rebuild)
      .def("connect", &VerletList::connect)
      .def("disconnect", &VerletList::disconnect)
    
      .def("getVerletCutoff", &VerletList::getVerletCutoff)
      ;
  }

}
