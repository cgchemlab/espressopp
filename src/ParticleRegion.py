#  Copyright (C) 2016
#      Jakub Krajniak (c) (jkrajniak at gmail.com)
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
****************************
**espressopp.ParticleRegion**
****************************


.. function:: espressopp.ParticleRegion(storage)

		:param storage: 
		:type storage:

.. function:: espressopp.ParticleRegion.has(pid)

		:param pid: 
		:type pid: 
		:rtype: 

.. function:: espressopp.ParticleRegion.show()

		:rtype: 

.. function:: espressopp.ParticleRegion.size()

		:rtype: 
"""
import _espressopp
import esutil
import pmi
from espressopp.esutil import cxxinit


class ParticleRegionLocal(_espressopp.ParticleRegion):
    def __init__(self, storage, left_bottom, right_top):
        if pmi.workerIsActive():
            cxxinit(self, _espressopp.ParticleRegion, storage)
            self.cxxclass.define_region(self, left_bottom, right_top)

    def add_type_id(self, type_id):
        if pmi.workerIsActive():
            self.cxxclass.add_type_id(self, type_id)

    def remove_type_id(self, type_id):
        if pmi.workerIsActive():
            self.cxxclass.remove_type_id(self, type_id)

    def show(self):
        if pmi.workerIsActive():
            self.cxxclass.show(self)

    def has(self, pid):
        if pmi.workerIsActive():
            return self.cxxclass.has(self, pid)

    def size(self):
        if pmi.workerIsActive():
            return self.cxxclass.size(self)


if pmi.isController:
    class ParticleRegion(object):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls='espressopp.ParticleRegionLocal',
            pmiinvoke=['get_particle_ids'],
            pmicall=['show', 'has', 'size', 'define_region', 'add_type_id', 'remove_type_id']
            )

