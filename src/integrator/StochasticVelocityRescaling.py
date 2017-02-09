#  Copyright (C) 2012,2013
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
*************************************************
espressopp.integrator.StochasticVelocityRescaling
*************************************************


.. function:: espressopp.integrator.StochasticVelocityRescaling(system)

		:param system: 
		:type system: 
"""
from espressopp.esutil import cxxinit
from espressopp import pmi

from espressopp.integrator.Extension import *
from _espressopp import integrator_StochasticVelocityRescaling 

class StochasticVelocityRescalingLocal(ExtensionLocal, integrator_StochasticVelocityRescaling):

    def __init__(self, system):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, integrator_StochasticVelocityRescaling, system)

    def add_valid_type_id(self, type_id):
        if pmi.workerIsActive():
            self.cxxclass.add_valid_type_id(self, type_id)

    def add_valid_types(self, types):
        if pmi.workerIsActive():
            for type_id in types:
                self.cxxclass.add_valid_type_id(self, type_id)

    def remove_valid_type_id(self, type_id):
        if pmi.workerIsActive():
            self.cxxclass.remove_valid_type_id(self, type_id)

if pmi.isController :
    class StochasticVelocityRescaling(Extension):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls =  'espressopp.integrator.StochasticVelocityRescalingLocal',
            pmiproperty = [ 'temperature', 'coupling' ],
            pmicall = ['add_valid_type_id', 'remove_valid_type_id', 'add_valid_types']
        )
