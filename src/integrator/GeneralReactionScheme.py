#  Copyright (C) 2018
#      Zidan Zhang (zidan.zhang@kuleuven.be)
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
***********************************************
**espressopp.integrator.GeneralReactionScheme**
***********************************************

This extension allows chainging properties of particles whenever they enters the defined geometrical region.

.. function:: espressopp.integrator.GeneralReactionScheme(system)

		:param system: 
		:type system:

"""
from espressopp.esutil import cxxinit
from espressopp import pmi
from espressopp.integrator.Extension import *
from _espressopp import integrator_GeneralReactionScheme


class GeneralReactionSchemeLocal(ExtensionLocal, integrator_GeneralReactionScheme):
    def __init__(self, system, interval):
        if pmi.workerIsActive():
            cxxinit(self, integrator_GeneralReactionScheme, system, interval)


if pmi.isController :
    class GeneralReactionScheme(Extension):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls =  'espressopp.integrator.GeneralReactionSchemeLocal',
        )
