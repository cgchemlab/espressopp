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
*****************************************
**espressopp.interaction.MultiTabulated**
*****************************************

This class of non-bonded tabulated potential allows to select table depends on the chemical conversion.
First step is to register table for certain range of conversion by method `register_table`

.. function:: espressopp.interaction.MultiTabulated(itype, cutoff)

		:param cutoff: (default: infinity)
		:type cutoff: float

.. function:: espressopp.interaction.MultiTabuleted.register_table(filename, itype, chem_conv_obs, min_value, max_value, default=False)

        :param filename: The filename with the table.
        :type filename: str
        :param itype: The type of interpolation, 1 - linear, 2 - Akima, 3 - Cubic
        :type itype: int
        :param chem_conv_obs: The instance of analysis.ChemicalConversion object.
        :type chem_conv_obs: espressopp.analysis.ChemicalConversion
        :param min_value: The minimum value of the conversion to use registered table.
        :type min_value: float
        :param max_value: The maximum value of the conversion to use registered table.
        :type max_value: float
        :param default: If set to true, this table will be used if for given conversion there is no other registered table
        :type default: bool

.. function:: espressopp.interaction.VerletListMultiTabulated(vl)

		:param vl: The VerletList object.
		:type vl: espressopp.VerletList.

"""
from espressopp import pmi, infinity
from espressopp.esutil import *

from espressopp.interaction.Potential import *
from espressopp.interaction.Interaction import *
from _espressopp import interaction_MultiTabulated, interaction_VerletListMultiTabulated


class MultiTabulatedLocal(PotentialLocal, interaction_MultiTabulated):
    def __init__(self, cutoff=infinity):
        if pmi.workerIsActive():
            cxxinit(self, interaction_MultiTabulated, cutoff)

    def register_table(self, filename, itype, chem_conv_obs, min_value, max_value, default=False):
        if pmi.workerIsActive():
            self.cxxclass.register_table(self, filename, itype, chem_conv_obs, min_value, max_value, default)

class VerletListMultiTabulatedLocal(InteractionLocal, interaction_VerletListMultiTabulated):
    def __init__(self, vl):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, interaction_VerletListMultiTabulated, vl)

    def setPotential(self, type1, type2, potential):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setPotential(self, type1, type2, potential)

    def getPotential(self, type1, type2):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            return self.cxxclass.getPotential(self, type1, type2)

if pmi.isController:
    class MultiTabulated(Potential):
        'The MultiTabulated potential.'
        pmiproxydefs = dict(
            cls='espressopp.interaction.MultiTabulatedLocal',
            pmicall=['register_table']
        )

    class VerletListMultiTabulated(Interaction):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls='espressopp.interaction.VerletListMultiTabulatedLocal',
            pmicall=['setPotential', 'getPotential']
        )