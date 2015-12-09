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
import pyh5md
import numpy as np


class DumpTopologyLocal(ParticleAccessLocal, io_DumpTopology):
    def __init__(self, system, integrator, h5md_file):
        if pmi.workerIsActive():
            cxxinit(self, io_DumpTopology, system, integrator)
            self.h5md_file = h5md_file
            self.tuple_index = 0
            self.tuple_data = {}
            if 'connectivity' not in self.h5md_file.file.f:
                self.h5md_file.file.f.create_group('connectivity')
            self.chunk_size = 256
            self.dt = integrator.dt

    def dump(self):
        if pmi.workerIsActive():
            self.cxxclass.dump(self)

    def observe_tuple(self, fpl, name, particle_group='atoms'):
        if pmi.workerIsActive():
            self.cxxclass.observe_tuple(self, fpl)
            g = pyh5md.TimeData(
                self.h5md_file.file.f['/connectivity'],
                name,
                shape=(None, None, 2),
                dtype=np.int,
                fillvalue=-1)
            g.attrs['particle_group'] = particle_group
            self.tuple_data[self.tuple_index] = g
            self.tuple_index += 1

    def update(self):
        if pmi.workerIsActive():
            raw_data = self.cxxclass.get_data(self)
            while raw_data != []:
                step = raw_data.pop()
                fpl_idx = raw_data.pop()
                fpl_size = raw_data.pop()
                if fpl_size == 0:
                    continue
                data = [
                    (raw_data.pop(), raw_data.pop())
                    for _ in range(fpl_size)
                ]
                g = self.tuple_data[fpl_idx]
                g.append(data, step, step*self.dt)

            self.cxxclass.clear_buffer(self)

if pmi.isController:
    class DumpTopology(ParticleAccess):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls =  'espressopp.io.DumpTopologyLocal',
            pmicall = ['dump', 'clear_buffer', 'observe_tuple', 'update'],
            pmiproperty = [],
            pmiinvoke = ['get_data']
        )
