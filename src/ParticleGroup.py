#  Copyright (C) 2012,2013,2015
#      Max Planck Institute for Polymer Research
#  Copyright (C) 2008,2009,2010,2011
#      Max-Planck-Institute for Polymer Research & Fraunhofer SCAI
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
************************
espressopp.ParticleGroup
************************


.. function:: espressopp.ParticleGroup(storage)

		:param storage: 
		:type storage: 

.. function:: espressopp.ParticleGroup.add(pid)

		:param pid: 
		:type pid: 
		:rtype: 

.. function:: espressopp.ParticleGroup.has(pid)

		:param pid: 
		:type pid: 
		:rtype: 

.. function:: espressopp.ParticleGroup.show()

		:rtype: 

.. function:: espressopp.ParticleGroup.size()

		:rtype: 
"""
import _espressopp
import esutil
import pmi
from espressopp.esutil import cxxinit

class ParticleGroupLocal(_espressopp.ParticleGroup):


    def __init__(self, storage):
        if pmi.workerIsActive():
            cxxinit(self, _espressopp.ParticleGroup, storage)

    def add(self, pid):
        if pmi.workerIsActive():
            self.cxxclass.add(self, pid)

    def show(self):
        if pmi.workerIsActive():
            self.cxxclass.show(self)

    def has(self, pid):
        if pmi.workerIsActive():
            return self.cxxclass.has(self, pid)

    def size(self):
        if pmi.workerIsActive():
            return self.cxxclass.size(self)


class ParticleGroupByTypeLocal(_espressopp.ParticleGroupByType):

    def __init__(self, storage, integrator):
        if pmi.workerIsActive():
            cxxinit(self, _espressopp.ParticleGroupByType, storage, integrator)

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
    class ParticleGroup(object):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls = 'espressopp.ParticleGroupLocal',
            pmicall = [ "add", "show", "has", "size" ]
            )

    class ParticleGroupByType(object):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls = 'espressopp.ParticleGroupByTypeLocal',
            pmiinvoke=['get_particle_ids'],
            pmicall = [ 'add_type_id', 'remove_type_id', 'show', 'has', 'size']
        )
