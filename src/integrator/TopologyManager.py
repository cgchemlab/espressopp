#
#  Copyright (C) 2015
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
*********************************
**espressopp.integrator.TopologyManager**
*********************************


.. function:: espressopp.integrator.TopologyManager(system)

		:param system: 
		:param verletlist: 
		:type system: 
		:type verletlist: 
"""
from espressopp.esutil import cxxinit
from espressopp import pmi

from _espressopp import integrator_TopologyManager

class TopologyManagerLocal(integrator_TopologyManager):

    def __init__(self, system):
        if pmi.workerIsActive():
            cxxinit(self, integrator_TopologyManager, system)

    def rebuild(self):
        if pmi.workerIsActive():
            self.cxxclass.rebuild(self)

    def observe_tuple(self, fpl):
        if pmi.workerIsActive():
            self.cxxclass.observe_tuple(self, fpl)

    def register_tuple(self, fpl, type1, type2=None):
        if pmi.workerIsActive():
            if type2 is None:
                type2 = type1
            self.cxxclass.register_tuple(self, fpl, type1, type2, 4)

    def register_triplet(self, ftl, type1, type2=None, type3=None):
        if pmi.workerIsActive():
            if type2 is None and type3 is None:
                type2 = type3 = type1
            self.cxxclass.register_triple(self, ftl, type1, type2, type3)

    def register_quadruplet(self, fql, type1, type2=None, type3=None, type4=None):
        if pmi.workerIsActive():
            if None in [type2, type3, type4]:
                type2 = type3 = type4 = type1
            self.cxxclass.register_quadruple(self, fql, type1, type2, type3, type4)

    def initialize_topology(self):
        if pmi.workerIsActive():
            self.cxxclass.initialize(self)

if pmi.isController :
    class TopologyManager(object):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls =  'espressopp.integrator.TopologyManagerLocal',
            pmicall = ['rebuild', 'observe_tuple', 'register_tuple', 'register_triplet',
                       'register_quadruplet', 'initialize_topology', 'exchange_data'],
            pmiinvoke = ['print_topology', 'get_neighbour_lists']
            )
