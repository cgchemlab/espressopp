#include "python.hpp"
#include "boost/bind.hpp"

#include "pairs/Set.hpp"
#include "pairs/Computer.hpp"

using namespace espresso::pairs;

void 
Set::foreach(Computer &computer) {
  computer.prepare();
  try {
    foreach(boost::bind(&Computer::apply, &computer, _1, _2, _3));
  } catch (ForeachBreak) {
    computer.finalize();
    throw;
  }
  computer.finalize();
}

void
Set::foreach(ConstComputer &computer) const {
  computer.prepare();
  try {
    foreach(boost::bind(&ConstComputer::apply, &computer, _1, _2, _3));
  } catch (ForeachBreak) {
    computer.finalize();
    throw;
  }
  computer.finalize();
}

void
Set::foreach(const Computer::SelfPtr computer) {
  foreach(*computer);
}

//////////////////////////////////////////////////
// REGISTRATION WITH PYTHON
//////////////////////////////////////////////////
void
Set::registerPython() {
  using namespace espresso::python;

  void (Set::*py_foreach)(const Computer::SelfPtr computer) 
    = &Set::foreach;

  // also register the abstract class Set to make virtual functions available
  // be careful: boost::noncopyable must be used for abstract classes with pure routines
  // no_init must be used as the abstract class Set has no constructor

  class_< Set, boost::noncopyable >("pairs_Set", no_init)
    .def("foreach", py_foreach);
  ;
}


