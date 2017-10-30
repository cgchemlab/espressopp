#  Copyright (C) 2017
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
********************************
**espressopp.analysis.BoxSize**
********************************

Returns the observable box size vector. Useful for NPT simulations.


.. function:: espressopp.analysis.BoxSize(system)

		:param system: 
		:type system: 
"""
from espressopp.esutil import cxxinit
from espressopp import pmi

from espressopp.analysis.Observable import *
from _espressopp import analysis_BoxSize

class BoxSizeLocal(ObservableLocal, analysis_BoxSize):
    def __init__(self, system):
        if pmi.workerIsActive():
            cxxinit(self, analysis_BoxSize, system)

if pmi.isController :
    class BoxSize(Observable):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls =  'espressopp.analysis.BoxSizeLocal'
            )
