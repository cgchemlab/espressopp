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
        """
        """
        if pmi.workerIsActive():
            self.cxxclass.rebuild(self)

    def observe(self, fpl):
        if pmi.workerIsActive():
            self.cxxclass.observe(self, fpl)

if pmi.isController :
    class TopologyManager(object):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls =  'espressopp.integrator.TopologyManagerLocal',
            pmicall = ['rebuild', 'observe']
            )
