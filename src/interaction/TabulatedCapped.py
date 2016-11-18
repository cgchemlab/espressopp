#  Copyright (C) 2016
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
**********************************************
**espressopp.interaction.TabulatedCapped**
**********************************************


.. function:: espressopp.interaction.TabulatedCapped(itype, filename, cutoff)

		:param itype: 
		:param filename: 
		:param cutoff: (default: infinity)
		:type itype: 
		:type filename: 
		:type cutoff: 

.. function:: espressopp.interaction.VerletListAdressTabulatedCapped(vl, fixedtupleList)

		:param vl: 
		:param fixedtupleList: 
		:type vl: 
		:type fixedtupleList: 

.. function:: espressopp.interaction.VerletListAdressTabulatedCapped.setPotentialAT(type1, type2, potential)

		:param type1: 
		:param type2: 
		:param potential: 
		:type type1: 
		:type type2: 
		:type potential: 

.. function:: espressopp.interaction.VerletListAdressTabulatedCapped.setPotentialCG(type1, type2, potential)

		:param type1: 
		:param type2: 
		:param potential: 
		:type type1: 
		:type type2: 
		:type potential: 

.. function:: espressopp.interaction.VerletListHadressTabulatedCapped(vl, fixedtupleList)

		:param vl: 
		:param fixedtupleList: 
		:type vl: 
		:type fixedtupleList: 

.. function:: espressopp.interaction.VerletListHadressTabulatedCapped.setPotentialAT(type1, type2, potential)

		:param type1: 
		:param type2: 
		:param potential: 
		:type type1: 
		:type type2: 
		:type potential: 

.. function:: espressopp.interaction.VerletListHadressTabulatedCapped.setPotentialCG(type1, type2, potential)

		:param type1: 
		:param type2: 
		:param potential: 
		:type type1: 
		:type type2: 
		:type potential: 

.. function:: espressopp.interaction.VerletListTabulatedCapped(vl)

		:param vl: 
		:type vl: 

.. function:: espressopp.interaction.VerletListTabulatedCapped.getPotential(type1, type2)

		:param type1: 
		:param type2: 
		:type type1: 
		:type type2: 
		:rtype: 

.. function:: espressopp.interaction.VerletListTabulatedCapped.setPotential(type1, type2, potential)

		:param type1: 
		:param type2: 
		:param potential: 
		:type type1: 
		:type type2: 
		:type potential: 

.. function:: espressopp.interaction.CellListTabulatedCapped(stor)

		:param stor: 
		:type stor: 

.. function:: espressopp.interaction.CellListTabulatedCapped.setPotential(type1, type2, potential)

		:param type1: 
		:param type2: 
		:param potential: 
		:type type1: 
		:type type2: 
		:type potential: 

.. function:: espressopp.interaction.FixedPairListTabulatedCapped(system, vl, potential)

		:param system: 
		:param vl: 
		:param potential: 
		:type system: 
		:type vl: 
		:type potential: 

.. function:: espressopp.interaction.FixedPairListTabulatedCapped.setPotential(potential)

		:param potential: 
		:type potential:

.. function:: espressopp.interaction.FixedPairListTypesTabulatedCapped(system, ftl)

        :param system: The Espresso++ system object.
        :type system: espressopp.System
        :param ftl: The FixedPair list.
        :type ftl: espressopp.FixedPairList

.. function:: espressopp.interaction.FixedPairListTypesTabulatedCapped.setPotential(type1, type2, potential)

        Defines bond potential for interaction between particles of types type1-type2-type3.

        :param type1: Type of particle 1.
        :type type1: int
        :param type2: Type of particle 2.
        :type type2: int
        :param potential: The potential to set up.
        :type potential: espressopp.interaction.Potential
"""
# -*- coding: iso-8859-1 -*-
from espressopp import pmi, infinity
from espressopp.esutil import *

from espressopp.interaction.Potential import *
from espressopp.interaction.Interaction import *
from _espressopp import interaction_TabulatedCapped, \
                      interaction_VerletListTabulatedCapped, \
                      interaction_VerletListAdressTabulatedCapped, \
                      interaction_VerletListHadressTabulatedCapped, \
                      interaction_VerletListDynamicResolutionTabulatedCapped, \
                      interaction_CellListTabulatedCapped, \
                      interaction_FixedPairListTabulatedCapped, \
                      interaction_FixedPairListTypesTabulatedCapped
from _espressopp import interaction_FixedPairListLambdaTabulatedCapped
from _espressopp import interaction_FixedPairListTypesLambdaTabulatedCapped

class TabulatedCappedLocal(PotentialLocal, interaction_TabulatedCapped):

    def __init__(self, itype, filename, cutoff=infinity, caprad=0.0):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, interaction_TabulatedCapped, int(itype), filename, cutoff, caprad)

class VerletListAdressTabulatedCappedLocal(InteractionLocal, interaction_VerletListAdressTabulatedCapped):

    def __init__(self, vl, fixedtupleList):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, interaction_VerletListAdressTabulatedCapped, vl, fixedtupleList)

    def setPotentialAT(self, type1, type2, potential):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setPotentialAT(self, type1, type2, potential)        
            
    def setPotentialCG(self, type1, type2, potential):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setPotentialCG(self, type1, type2, potential)
            
class VerletListHadressTabulatedCappedLocal(InteractionLocal, interaction_VerletListHadressTabulatedCapped):

    def __init__(self, vl, fixedtupleList):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, interaction_VerletListHadressTabulatedCapped, vl, fixedtupleList)

    def setPotentialAT(self, type1, type2, potential):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setPotentialAT(self, type1, type2, potential)        
            
    def setPotentialCG(self, type1, type2, potential):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setPotentialCG(self, type1, type2, potential)

class VerletListTabulatedCappedLocal(InteractionLocal, interaction_VerletListTabulatedCapped):

    def __init__(self, vl):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, interaction_VerletListTabulatedCapped, vl)

    def setPotential(self, type1, type2, potential):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setPotential(self, type1, type2, potential)

    def getPotential(self, type1, type2):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            return self.cxxclass.getPotential(self, type1, type2)

class VerletListDynamicResolutionTabulatedCappedLocal(InteractionLocal, interaction_VerletListDynamicResolutionTabulatedCapped):
    def __init__(self, vl, cg_form):
        if pmi.workerIsActive():
            cxxinit(self, interaction_VerletListDynamicResolutionTabulatedCapped, vl, cg_form)

    def setPotential(self, type1, type2, potential):
        if pmi.workerIsActive():
            self.cxxclass.setPotential(self, type1, type2, potential)

    def getPotential(self, type1, type2):
        if pmi.workerIsActive():
            return self.cxxclass.getPotential(self, type1, type2)

    def getVerletListLocal(self):
        if pmi.workerIsActive():
            return self.cxxclass.getVerletList(self)

    def setMaxForce(self, max_force):
        if pmi.workerIsActive():
            self.cxxclass.setMaxForce(self, max_force)

class CellListTabulatedCappedLocal(InteractionLocal, interaction_CellListTabulatedCapped):

    def __init__(self, stor):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, interaction_CellListTabulatedCapped, stor)

    def setPotential(self, type1, type2, potential):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setPotential(self, type1, type2, potential)


class FixedPairListTabulatedCappedLocal(InteractionLocal, interaction_FixedPairListTabulatedCapped):

    def __init__(self, system, vl, potential):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, interaction_FixedPairListTabulatedCapped, system, vl, potential)

    def setPotential(self, potential):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setPotential(self, potential)


class FixedPairListTypesTabulatedCappedLocal(InteractionLocal, interaction_FixedPairListTypesTabulatedCapped):
    def __init__(self, system, vl):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, interaction_FixedPairListTypesTabulatedCapped, system, vl)

    def setPotential(self, type1, type2, potential):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setPotential(self, type1, type2, potential)

    def getPotential(self, type1, type2):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            return self.cxxclass.getPotential(self, type1, type2)

    def setFixedPairList(self, fixedpairlist):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setFixedPairList(self, fixedpairlist)

    def getFixedPairList(self):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            return self.cxxclass.getFixedPairList(self)


class FixedPairListLambdaTabulatedCappedLocal(InteractionLocal, interaction_FixedPairListLambdaTabulatedCapped):

    def __init__(self, system, vl, potential):
        if pmi.workerIsActive():
            cxxinit(self, interaction_FixedPairListLambdaTabulatedCapped, system, vl, potential)

    def setPotential(self, potential):
        if pmi.workerIsActive():
            self.cxxclass.setPotential(self, potential)


class FixedPairListTypesLambdaTabulatedCappedLocal(InteractionLocal, interaction_FixedPairListTypesLambdaTabulatedCapped):
    def __init__(self, system, fpl):
        if pmi.workerIsActive():
            cxxinit(self, interaction_FixedPairListTypesLambdaTabulatedCapped, system, fpl)

    def setPotential(self, type1, type2, potential):
        if pmi.workerIsActive():
            self.cxxclass.setPotential(self, type1, type2, potential)

    def getPotential(self, type1, type2):
        if pmi.workerIsActive():
            return self.cxxclass.getPotential(self, type1, type2)

    def setFixedPairList(self, fixedpairlist):
        if pmi.workerIsActive():
            self.cxxclass.setFixedPairList(self, fixedpairlist)

    def getFixedPairList(self):
        if pmi.workerIsActive():
            return self.cxxclass.getFixedPairList(self)


if pmi.isController:
    class TabulatedCapped(Potential):
        'The TabulatedCapped potential.'
        pmiproxydefs = dict(
            cls = 'espressopp.interaction.TabulatedCappedLocal',
            pmiproperty = ['itype', 'filename', 'cutoff', 'caprad']
            )
        
    class VerletListAdressTabulatedCapped(Interaction):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls =  'espressopp.interaction.VerletListAdressTabulatedCappedLocal',
            pmicall = ['setPotentialAT', 'setPotentialCG']
            )
            
    class VerletListHadressTabulatedCapped(Interaction):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls =  'espressopp.interaction.VerletListHadressTabulatedCappedLocal',
            pmicall = ['setPotentialAT', 'setPotentialCG']
            )
        
    class VerletListTabulatedCapped(Interaction):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls =  'espressopp.interaction.VerletListTabulatedCappedLocal',
            pmicall = ['setPotential','getPotential']
            )

    class VerletListDynamicResolutionTabulatedCapped(Interaction):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls =  'espressopp.interaction.VerletListDynamicResolutionTabulatedCappedLocal',
            pmicall = ['setPotential', 'getPotential', 'getVerletList', 'setMaxForce']
        )

    class CellListTabulatedCapped(Interaction):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls =  'espressopp.interaction.CellListTabulatedCappedLocal',
            pmicall = ['setPotential']
            )
        
    class FixedPairListTabulatedCapped(Interaction):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls =  'espressopp.interaction.FixedPairListTabulatedCappedLocal',
            pmicall = ['setPotential', 'setFixedPairList', 'getFixedPairList']
            )

    class FixedPairListTypesTabulatedCapped(Interaction):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls =  'espressopp.interaction.FixedPairListTypesTabulatedCappedLocal',
            pmicall = ['setPotential','getPotential','setFixedPairList','getFixedPairList']
        )

    class FixedPairListLambdaTabulatedCapped(Interaction):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls =  'espressopp.interaction.FixedPairListLambdaTabulatedCappedLocal',
            pmicall = ['setPotential', 'setFixedPairList', 'getFixedPairList']
        )

    class FixedPairListTypesLambdaTabulatedCapped(Interaction):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls =  'espressopp.interaction.FixedPairListTypesLambdaTabulatedCappedLocal',
            pmicall = ['setPotential','getPotential','setFixedPairList','getFixedPairList']
        )