#  Copyright (C) 2012-2016
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
*****************************************
**FENE interaction**
*****************************************

Implementation of the Finitely Extensible Non-linear Elastic potential:

.. math:: 

        U(r) = -\frac{1}{2}r_{\mathrm{max}}^2  K \log\left[1 - \left(\frac{r - r_{0}}{r_{\mathrm{max}}}\right)^2\right]


.. function:: espressopp.interaction.FENE(K, r0, rMax, cutoff, shift)

                :param real K: (default: 1.0)
                :param real r0: (default: 0.0)
		:param rMax: (default: 1.0)
		:param cutoff: (default: infinity)
		:param shift: (default: 0.0)
		:type rMax: real
                :type cutoff: real
		:type shift: real

.. function:: espressopp.interaction.FixedPairListFENE(system, pair_list, potential)

                :param object system: your system :func:`espressopp.System`
                :param object pair_list: list of bonds  :func:`espressopp.FixedPairList`
                :param object potential: :func:`espressopp.interaction.FENE`

.. function:: espressopp.interaction.FixedPairListFENE.getFixedPairList()

		:rtype: A Python list of lists.

.. function:: espressopp.interaction.FixedPairListFENE.getPotential()

                :rtype: object

.. function:: espressopp.interaction.FixedPairListFENE.setFixedPairList(pair_list)

                :param pair_list:
                :type pair_list: fixedpairlist

.. function:: espressopp.interaction.FixedPairListFENE.setPotential(potential)

		:param potential: 
		:type potential: 

**Example of usage**

>>> # The following example shows how to bond particle 1 to particles 0 and 2 by a FENE potential.
>>> # We assume the particles are already in the storage of the system
>>> # Initialize list of pairs that will be bonded by FENE
>>> pair_list = espressopp.FixedPairList(system.storage)
>>> # Set which pairs belong to the pair_list i.e. particle 0 is bonded to particles 1 and 2.
>>> pair_list.addBonds([(0,1),(1,2)])
>>> # Initialize the potential and set up the parameters.
>>> potFENE   = espressopp.interaction.FENE(K=30.0, r0=0.0, rMax=1.5)
>>> # Set which system, pair list and potential is the interaction associated with.
>>> interFENE = espressopp.interaction.FixedPairListFENE(system, pair_list, potFENE)
>>> # Add the interaction to the system.
>>> system.addInteraction(interFENE)

"""

from espressopp import pmi, infinity
from espressopp.esutil import *

from espressopp.interaction.Potential import *
from espressopp.interaction.Interaction import *
from _espressopp import interaction_FENE, interaction_FixedPairListFENE
from _espressopp import interaction_FixedPairListTypesFENE
from _espressopp import interaction_FixedPairListLambdaFENE
from _espressopp import interaction_FixedPairListTypesLambdaFENE

class FENELocal(PotentialLocal, interaction_FENE):

    def __init__(self, K=1.0, r0=0.0, rMax=1.0, 
                 cutoff=infinity, shift=0.0):
        if pmi.workerIsActive():
            if shift == "auto":
                cxxinit(self, interaction_FENE, K, r0, rMax, cutoff)
            else:
                cxxinit(self, interaction_FENE, K, r0, rMax, cutoff, shift)

class FixedPairListFENELocal(InteractionLocal, interaction_FixedPairListFENE):

    def __init__(self, system, vl, potential):
        if pmi.workerIsActive():
            cxxinit(self, interaction_FixedPairListFENE, system, vl, potential)

    def setPotential(self, potential):
        if pmi.workerIsActive():
            self.cxxclass.setPotential(self, potential)

    def getPotential(self):
        if pmi.workerIsActive():
            return self.cxxclass.getPotential(self)

    def setFixedPairList(self, fixedpairlist):
        if pmi.workerIsActive():
            self.cxxclass.setFixedPairList(self, fixedpairlist)

    def getFixedPairList(self):
        if pmi.workerIsActive():
            return self.cxxclass.getFixedPairList(self)

class FixedPairListTypesFENELocal(InteractionLocal, interaction_FixedPairListTypesFENE):
    def __init__(self, system, vl):
        if pmi.workerIsActive():
            cxxinit(self, interaction_FixedPairListTypesFENE, system, vl)

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

class FixedPairListLambdaFENELocal(InteractionLocal, interaction_FixedPairListLambdaFENE):

    def __init__(self, system, fpl, potential):
        if pmi.workerIsActive():
            cxxinit(self, interaction_FixedPairListLambdaFENE, system, fpl, potential)

    def setPotential(self, potential):
        if pmi.workerIsActive():
            self.cxxclass.setPotential(self, potential)

    def setFixedPairList(self, fixedpairlist):
        if pmi.workerIsActive():
            self.cxxclass.setFixedPairList(self, fixedpairlist)


    def getFixedPairList(self):
        if pmi.workerIsActive():
            return self.cxxclass.getFixedPairList(self)

class FixedPairListTypesLambdaFENELocal(InteractionLocal, interaction_FixedPairListTypesLambdaFENE):
    def __init__(self, system, fpl):
        if pmi.workerIsActive():
            cxxinit(self, interaction_FixedPairListTypesLambdaFENE, system, fpl)

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
    class FENE(Potential):
        pmiproxydefs = dict(
            cls = 'espressopp.interaction.FENELocal',
            pmiproperty = ['K', 'r0', 'rMax']
            )

    class FixedPairListFENE(Interaction):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls =  'espressopp.interaction.FixedPairListFENELocal',
            pmicall = ['setPotential','getPotential','setFixedPairList', 'getFixedPairList']
            )

    class FixedPairListLambdaFENE(Interaction):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls='espressopp.interaction.FixedPairListLambdaFENELocal',
            pmicall = ['setPotential','getPotential','setFixedPairList','getFixedPairList']
        )

    class FixedPairListTypesFENE(Interaction):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls =  'espressopp.interaction.FixedPairListTypesFENELocal',
            pmicall = ['setPotential','getPotential','setFixedPairList','getFixedPairList']
        )

    class FixedPairListTypesLambdaFENE(Interaction):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls =  'espressopp.interaction.FixedPairListTypesLambdaFENELocal',
            pmicall = ['setPotential','getPotential','setFixedPairList','getFixedPairList']
        )