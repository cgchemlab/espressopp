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

#ifndef _IO_ReadNumpy_HPP
#define _IO_ReadNumpy_HPP

#include <algorithm>
#include <string>
#include "boost/python/numeric.hpp"

#include "esutil/Error.hpp"

#include "storage/Storage.hpp"
#include "bc/BC.hpp"

namespace espressopp {
namespace io {

class ReadNumpy: public SystemAccess {
 public:
  explicit ReadNumpy(shared_ptr <System> system) : SystemAccess(system) {
  }
  ~ReadNumpy() {
  }

  static void registerPython();
 protected:
  static LOG4ESPP_DECL_LOGGER(logger);
 private:
  long* LoadIds(PyObject *ids, long &num_particles);  // NOLINT
  bool LoadPosition(PyObject *ids, PyObject *data);

  Py_buffer LoadBuffer(PyObject *data);
};
}  // end namespace io
}  // end namespace espressopp

#endif
