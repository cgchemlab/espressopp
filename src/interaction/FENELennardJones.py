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
***************************************
espressopp.interaction.FENELennardJones
***************************************

Implementation of the Finitely Extensible Non-linear Elastic potential combined with Lennard-Jones potential.
It is used for bead-spring polymer model.

.. math:: 

        U(r) = -\frac{1}{2}r_{\mathrm{max}}^2  K \log\left[1 - \left(\frac{r - r_{0}}{r_{\mathrm{max}}}\right)^2\right]


.. function:: espressopp.interaction.FENELennardJones(K, r0, rMax, sigma, epsilon, cutoff, shift)

        :param real K: (default: 1.0)
        :param real r0: (default: 0.0)
		:param rMax: (default: 1.0)
		:param sigma: (default: 1.0)
		:param epsilon: (default 1.0)
		:param cutoff: (default: infinity)
		:param shift: (default: 0.0)
		:type K: real
		:type r0: real
		:type sigma: real
		:type epsilon: real
		:type rMax: real
        :type cutoff: real
		:type shift: real

.. function:: espressopp.interaction.FixedPairListFENELennardJones(system, pair_list, potential)

                :param object system: your system :func:`espressopp.System`
                :param object pair_list: list of bonds  :func:`espressopp.FixedPairList`
                :param object potential: :func:`espressopp.interaction.FENELennardJones`

.. function:: espressopp.interaction.FixedPairListFENELennardJones.getFixedPairList()

		:rtype: A Python list of lists.

.. function:: espressopp.interaction.FixedPairListFENELennardJones.getPotential()

                :rtype: object

.. function:: espressopp.interaction.FixedPairListFENELennardJones.setFixedPairList(pair_list)

                :param pair_list:
                :type pair_list: fixedpairlist

.. function:: espressopp.interaction.FixedPairListFENELennardJones.setPotential(potential)

		:param potential: 
		:type potential: 

**Example of usage**

>>> # The following example shows how to bond particle 1 to particles 0 and 2 by a FENELennardJones potential.
>>> # We assume the particles are already in the storage of the system
>>> # Initialize list of pairs that will be bonded by FENELennardJones
>>> pair_list = espressopp.FixedPairList(system.storage)
>>> # Set which pairs belong to the pair_list i.e. particle 0 is bonded to particles 1 and 2.
>>> pair_list.addBonds([(0,1),(1,2)])
>>> # Initialize the potential and set up the parameters.
>>> potFENELennardJones   = espressopp.interaction.FENELennardJones(K=30.0, r0=0.0, rMax=1.5)
>>> # Set which system, pair list and potential is the interaction associated with.
>>> interFENELennardJones = espressopp.interaction.FixedPairListFENELennardJones(system, pair_list, potFENELennardJones)
>>> # Add the interaction to the system.
>>> system.addInteraction(interFENELennardJones)

"""

from espressopp import pmi, infinity
from espressopp.esutil import *

from espressopp.interaction.Potential import *
from espressopp.interaction.Interaction import *
from _espressopp import interaction_FENELennardJones, interaction_FixedPairListFENELennardJones
from _espressopp import interaction_FixedPairListTypesFENELennardJones
from _espressopp import interaction_FixedPairListLambdaFENELennardJones
from _espressopp import interaction_FixedPairListTypesLambdaFENELennardJones


class FENELennardJonesLocal(PotentialLocal, interaction_FENELennardJones):

    def __init__(self, K=1.0, r0=0.0, rMax=1.0,
                 sigma=1.0, epsilon=1.0, cutoff=infinity, shift=0.0):
        if pmi.workerIsActive():
            if shift == "auto":
                cxxinit(self, interaction_FENELennardJones,
                        K, r0, rMax, sigma, cutoff)
            else:
                cxxinit(self, interaction_FENELennardJones, K,
                        r0, rMax, sigma, epsilon, cutoff, shift)


class FixedPairListFENELennardJonesLocal(InteractionLocal, interaction_FixedPairListFENELennardJones):

    def __init__(self, system, vl, potential):
        if pmi.workerIsActive():
            cxxinit(self, interaction_FixedPairListFENELennardJones,
                    system, vl, potential)

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


class FixedPairListTypesFENELennardJonesLocal(InteractionLocal, interaction_FixedPairListTypesFENELennardJones):
    def __init__(self, system, vl):
        if pmi.workerIsActive():
            cxxinit(self, interaction_FixedPairListTypesFENELennardJones, system, vl)

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


class FixedPairListLambdaFENELennardJonesLocal(InteractionLocal, interaction_FixedPairListLambdaFENELennardJones):

    def __init__(self, system, fpl, potential):
        if pmi.workerIsActive():
            cxxinit(self, interaction_FixedPairListLambdaFENELennardJones,
                    system, fpl, potential)

    def setPotential(self, potential):
        if pmi.workerIsActive():
            self.cxxclass.setPotential(self, potential)

    def setFixedPairList(self, fixedpairlist):
        if pmi.workerIsActive():
            self.cxxclass.setFixedPairList(self, fixedpairlist)

    def getFixedPairList(self):
        if pmi.workerIsActive():
            return self.cxxclass.getFixedPairList(self)


class FixedPairListTypesLambdaFENELennardJonesLocal(InteractionLocal, interaction_FixedPairListTypesLambdaFENELennardJones):
    def __init__(self, system, fpl):
        if pmi.workerIsActive():
            cxxinit(
                self, interaction_FixedPairListTypesLambdaFENELennardJones, system, fpl)

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
    class FENELennardJones(Potential):
        pmiproxydefs = dict(
            cls='espressopp.interaction.FENELennardJonesLocal',
            pmiproperty=['K', 'r0', 'rMax', 'sigma', 'epsilon']
        )

    class FixedPairListFENELennardJones(Interaction):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls='espressopp.interaction.FixedPairListFENELennardJonesLocal',
            pmicall=['setPotential', 'getPotential',
                     'setFixedPairList', 'getFixedPairList']
        )

    class FixedPairListLambdaFENELennardJones(Interaction):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls='espressopp.interaction.FixedPairListLambdaFENELennardJonesLocal',
            pmicall=['setPotential', 'getPotential',
                     'setFixedPairList', 'getFixedPairList']
        )

    class FixedPairListTypesFENELennardJones(Interaction):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls='espressopp.interaction.FixedPairListTypesFENELennardJonesLocal',
            pmicall=['setPotential', 'getPotential',
                     'setFixedPairList', 'getFixedPairList']
        )

    class FixedPairListTypesLambdaFENELennardJones(Interaction):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls='espressopp.interaction.FixedPairListTypesLambdaFENELennardJonesLocal',
            pmicall=['setPotential', 'getPotential',
                     'setFixedPairList', 'getFixedPairList']
        )
