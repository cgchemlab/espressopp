#  Copyright (C) 2016
#      Jakub Krajniak (jkrajniak at gmail.com)
#
#  This file is part of ESPResSo++.
#
#  ESPResSo++ is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  ESPResSo++ is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program.  If not, see <http://www.gnu.org/licenses/>.


r"""
*********************************************
**ReadNumpy** - IO Object
*********************************************
"""

from _espressopp import io_ReadNumpy
from espressopp import pmi
from espressopp.esutil import cxxinit


class ReadNumpyLocal(io_ReadNumpy):
    def __init__(self, system):
        cxxinit(self, io_ReadNumpy, system)

    def load_position(self, ids, input_ndarray):
        if pmi.workerIsActive():
            self.cxxclass.load_position(self, ids, input_ndarray)


if pmi.isController:
    class ReadNumpy():
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls='espressopp.io.ReadNumpyLocal',
            pmicall=['load_position'],
        )
