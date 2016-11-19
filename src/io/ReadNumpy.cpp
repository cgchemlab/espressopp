/*
  Copyright (C) 2016
      Jakub Krajniak (jkrajniak at gmail.com)

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

#include <Python.h>

#include "ReadNumpy.hpp"
#include "python.hpp"


using namespace espressopp;  // NOLINT

namespace espressopp {
namespace io {


Py_buffer ReadNumpy::LoadBuffer(PyObject *data) {
  Py_buffer view_data;

  if (PyObject_CheckBuffer(data) == 0) {
    std::stringstream msg;
    msg << "Data does not have py_buffer interface";
    throw std::runtime_error(msg.str());
  }

  if (PyObject_GetBuffer(data, &view_data, PyBUF_ANY_CONTIGUOUS | PyBUF_FORMAT) == -1) {
    throw std::runtime_error("Wrong format of data");
  }

  return view_data;
}

long* ReadNumpy::LoadIds(PyObject *ids, long &num_particles) {  // NOLINT
  Py_buffer view_data = LoadBuffer(ids);
  if (strcmp(view_data.format, "l") != 0) {
    throw std::runtime_error("Expected an array of integers");
  }
  if (view_data.ndim != 1) {
    throw std::runtime_error("Expected a 1-d array");
  }

  long *buf = static_cast<long*>(view_data.buf);  //  NOLINT
  num_particles = view_data.len / view_data.itemsize;
  long *particle_ids = new long[num_particles];  // NOLINT
  std::copy(buf, buf+num_particles, particle_ids);

  return particle_ids;
}

bool ReadNumpy::LoadPosition(PyObject *ids, PyObject *data) {
  Py_buffer view_data = LoadBuffer(data);

  long *particle_ids;  // NOLINT
  long num_particles = 0;  // NOLINT
  particle_ids = LoadIds(ids, num_particles);

  if (strcmp(view_data.format, "d") != 0) {
    throw std::runtime_error("Expected an array of doubles");
  }

  if (view_data.ndim != 2) {
    throw std::runtime_error("Expected a 2-d array");
  }

  if ((view_data.len / view_data.itemsize) != (3*num_particles)) {
    throw std::runtime_error("Wrong number of elements");
  }

  System &system = getSystemRef();
  shared_ptr<storage::Storage> storage = system.storage;

  double *buf = static_cast<double*>(view_data.buf);
  for (longint i = 0; i < num_particles; i++) {
    long pid = particle_ids[i];  //NOLINT

    Particle *p = system.storage->lookupRealParticle(pid);
    if (p) {
      Real3D &pos = p->position();
      for (int j = 0; j < 3; j++) {
        pos[j] = buf[i*3+j];
      }
    }
  }

  delete[] particle_ids;
  return true;
}

// Python wrapping
void ReadNumpy::registerPython() {
  using namespace espressopp::python;  // NOLINT

  boost::python::numeric::array::set_module_and_type("numpy", "ndarray");

  class_< ReadNumpy>
      ("io_ReadNumpy", init<shared_ptr<System> >())
          .def("load_position", &ReadNumpy::LoadPosition);
}

}  // end namespace io
}  // end namespace espressopp
