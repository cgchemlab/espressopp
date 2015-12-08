#  Copyright (c) 2015
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
**DumpPairs** - IO Object
*********************************************

* `dump()`

  Properties

* `h5md_file`
  HDF5 file object.
"""

from espressopp.esutil import cxxinit
from espressopp import pmi

from espressopp.ParticleAccess import *
from _espressopp import io_DumpTopology

import collections


class DumpTopologyLocal(ParticleAccessLocal, io_DumpTopology):
    def __init__(self, system, integrator):
        if pmi.workerIsActive():
            cxxinit(self, io_DumpTopology, system, integrator)
            self.tuple_index = 0
            self.tuple_data = collections.defaultdict(dict)
  
    def dump(self):
        if pmi.workerIsActive():
            self.cxxclass.dump(self)

    def observe_tuple(self, fpl, name, particle_group='atoms'):
        if pmi.workerIsActive():
            selx.cxxclass.observe_tuple(self, fpl)
            self.tuple_data[particle_group][name] = fpl
            self.tuple_index += 1

    def update(self):
        pass
  
  
if pmi.isController :
    class DumpTopology(ParticleAccess):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls =  'espressopp.io.DumpTopologyLocal',
            pmicall = ['dump', 'clear_buffer', 'observe_tuple'],
            pmiproperty = [],
            pmiinvoke = ['get_data']
        )