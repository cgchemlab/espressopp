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
************************************************************
**espressopp.interaction.AngularHarmonic**
************************************************************

Calculates the Angular Harmonic interaction

.. math::

	U = K (\theta - \theta_0)^2

.. function:: espressopp.interaction.AngularHarmonic(K, theta0)

		:param K: (default: 1.0)
		:param theta0: (default: 0.0)
		:type K: real
		:type theta0: real

.. function:: espressopp.interaction.FixedTripleListAngularHarmonic(system, vl, potential)

		:param system: 
		:param vl: 
		:param potential: 
		:type system: 
		:type vl: 
		:type potential: 

.. function:: espressopp.interaction.FixedTripleListAngularHarmonic.setPotential(type1, type2, potential)

		:param type1: 
		:param type2: 
		:param potential: 
		:type type1: 
		:type type2: 
		:type potential: 
"""
from espressopp import pmi
from espressopp.esutil import *

from espressopp.interaction.AngularPotential import *
from espressopp.interaction.Interaction import *
from _espressopp import interaction_AngularHarmonic, interaction_FixedTripleListAngularHarmonic
from _espressopp import interaction_FixedTripleListTypesAngularHarmonic

class AngularHarmonicLocal(AngularPotentialLocal, interaction_AngularHarmonic):

    def __init__(self, K=1.0, theta0=0.0):

        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, interaction_AngularHarmonic, K, theta0)

class FixedTripleListAngularHarmonicLocal(InteractionLocal, interaction_FixedTripleListAngularHarmonic):

    def __init__(self, system, vl, potential):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, interaction_FixedTripleListAngularHarmonic, system, vl, potential)

    def setPotential(self, type1, type2, potential):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setPotential(self, type1, type2, potential)

class FixedTripleListTypesAngularHarmonicLocal(InteractionLocal, interaction_FixedTripleListTypesAngularHarmonic):
    def __init__(self, system, vl):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, interaction_FixedTripleListTypesAngularHarmonic, system, vl)

    def setPotential(self, type1, type2, type3, potential):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setPotential(self, type1, type2, type3, potential)

    def getPotential(self, type1, type2):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            return self.cxxclass.getPotential(self, type1, type2)

    def setFixedTripleList(self, fixedtriplelist):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setFixedTripleList(self, fixedtriplelist)

    def getFixedTripleList(self):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            return self.cxxclass.getFixedTripleList(self)

if pmi.isController:
    class AngularHarmonic(AngularPotential):
        'The AngularHarmonic potential.'
        pmiproxydefs = dict(
            cls = 'espressopp.interaction.AngularHarmonicLocal',
            pmiproperty = ['K', 'theta0']
            )

    class FixedTripleListAngularHarmonic(Interaction):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls =  'espressopp.interaction.FixedTripleListAngularHarmonicLocal',
            pmicall = ['setPotential', 'getFixedTripleList']
            )
    class FixedTripleListTypesAngularHarmonic(Interaction):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls =  'espressopp.interaction.FixedTripleListTypesAngularHarmonicLocal',
            pmicall = ['setPotential','getPotential','setFixedTripleList','getFixedTripleList']
        )