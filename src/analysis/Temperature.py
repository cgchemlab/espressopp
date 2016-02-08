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
***********************************
**espressopp.analysis.Temperature**
***********************************


.. function:: espressopp.analysis.Temperature(system)

		:param system: 
		:type system: 
"""
from espressopp.esutil import cxxinit
from espressopp import pmi

from espressopp.analysis.Observable import *
from _espressopp import analysis_Temperature

class TemperatureLocal(ObservableLocal, analysis_Temperature):

    def __init__(self, system):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, analysis_Temperature, system)
    
    def add_type(self, type_id):
        if pmi.workerIsActive():
            self.cxxclass.add_type(self, type_id)

    def remove_type(self, type_id):
        if pmi.workerIsActive():
            ret_val = self.cxxclass.remove_type(self, type_id)
            if not ret_val:
                print('Warning, type {} not found'.format(type_id))


if pmi.isController :
    class Temperature(Observable):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls =  'espressopp.analysis.TemperatureLocal',
            pmicall = ['add_type', 'remove_type']
            )
