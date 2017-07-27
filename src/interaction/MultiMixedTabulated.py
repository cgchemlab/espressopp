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
**********************************************
**espressopp.interaction.MultiMixedTabulated**
**********************************************

This class of non-bonded tabulated potential allows to select table depends on the chemical conversion.
First step is to register table for certain range of conversion by method `register_table`

This kind of potential effectively mix two tabulated potentials, based
on the input conversion value :math:`x`

.. math::
    U_{eff} = xU_I + (1-x)U_{II}

where :math:`x` is a parameter that can change during the simulation. :math:`U_I` is table_1 and :math:`U_{II}` is table_2.

For given conversion :math:`val` the table1, table2 are selected from the registered tables and are used to calculated effective force and
energy.

The table1,table2 is selected iff.

.. math::
    p2 > val <= p2

where :math:`val` is the current conversion and :math:`p1`,:math:`p2` are the set condtions.

.. function:: espressopp.interaction.MultiMixedTabulated(itype, cutoff)

        :param itype: The type of interpolation, 1 - linear, 2 - Akima, 3 - Cubic
        :type itype: int
		:param cutoff: (default: infinity)
		:type cutoff: float

.. function:: espressopp.interaction.MultiTabulated.register_table(filename, chem_conv_obs, min_value, max_value)

        :param filename: The filename with the table.
        :type filename: str
        :param chem_conv_obs: The instance of analysis.ChemicalConversion object.
        :type chem_conv_obs: espressopp.analysis.ChemicalConversion
        :param min_value: The minimum value of the conversion to use registered table.
        :type min_value: float
        :param max_value: The maximum value of the conversion to use registered table.
        :type max_value: float

.. function:: espressopp.interaction.VerletListMultiMixedTabulated(vl)

		:param vl: The VerletList object.
		:type vl: espressopp.VerletList.

"""
from espressopp import pmi, infinity
from espressopp.esutil import *

from espressopp.interaction.Potential import *
from espressopp.interaction.Interaction import *
from _espressopp import interaction_MultiMixedTabulated
from _espressopp import interaction_VerletListMultiMixedTabulated
from _espressopp import interaction_VerletListDynamicResolutionMultiMixedTabulated


class MultiMixedTabulatedLocal(PotentialLocal, interaction_MultiMixedTabulated):
    def __init__(self, itype, cutoff=infinity):
        if pmi.workerIsActive():
            cxxinit(self, interaction_MultiMixedTabulated, itype, cutoff)

    def register_table(self, tab1, tab2, chem_conv_obs, min_value, max_value):
        if pmi.workerIsActive():
            self.cxxclass.register_table(self, tab1, tab2, chem_conv_obs, min_value, max_value)

class VerletListMultiMixedTabulatedLocal(InteractionLocal, interaction_VerletListMultiMixedTabulated):
    def __init__(self, vl):
        if pmi.workerIsActive():
            cxxinit(self, interaction_VerletListMultiMixedTabulated, vl)

    def setPotential(self, type1, type2, potential):
        if pmi.workerIsActive():
            self.cxxclass.setPotential(self, type1, type2, potential)

    def getPotential(self, type1, type2):
        if pmi.workerIsActive():
            return self.cxxclass.getPotential(self, type1, type2)

if pmi.isController:
    class MultiMixedTabulated(Potential):
        pmiproxydefs = dict(
            cls='espressopp.interaction.MultiMixedTabulatedLocal',
            pmicall=['register_table']
        )

    class VerletListMultiMixedTabulated(Interaction):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls='espressopp.interaction.VerletListMultiMixedTabulatedLocal',
            pmicall=['setPotential', 'getPotential']
        )