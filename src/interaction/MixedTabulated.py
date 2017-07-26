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
**Mixed tabulated non-bonded potentials**
*****************************************

This kind of potential effectively mix two tabulated potentials, based
on the input parameter.

.. math::
    U_{eff} = xU_I + (1-x)U_{II}

where :math:`x` is a parameter that can change during the simulation. U_I is table_1 and U_II is table_2.


.. function:: espressopp.interaction.MixedTabulated(itype, table1, table2, conversion, cutoff)

        :param itype: The interpolation type (1 - linear, 2 - spline, 3 - cubic)
        :type itype: int
        :param table1: The filename of tabulated potential I
        :type table1: str
        :param tabl2: The filename of tabulated potential II
        :type table2: str
        :param conversion: The espressopp.analysis.ChemicalConversion object.
        :type conversion: espressopp.analysis.ChemicalConversion
        :param mix_value: The default initial :math:`x` value.
        :type mix_value: float
        :param cutoff: (default: infinity)
        :type cutoff: float

.. function:: espressopp.interaction.VerletListMixedTabulated(vl)

		:param vl: The VerletList object.
		:type vl: espressopp.VerletList.

"""
from espressopp import pmi, infinity
from espressopp.esutil import *

from espressopp.interaction.Potential import *
from espressopp.interaction.Interaction import *
from _espressopp import interaction_MixedTabulated, interaction_VerletListMixedTabulated


class MixedTabulatedLocal(PotentialLocal, interaction_MixedTabulated):
    def __init__(self, itype, table1, table2, conversion=None, mix_value=1.0, cutoff=infinity):
        if pmi.workerIsActive():
            if conversion:
                cxxinit(self, interaction_MixedTabulated, itype, table1, table2, conversion, mix_value, cutoff)
            else:
                cxxinit(self, interaction_MixedTabulated, itype, table1, table2, mix_value, cutoff)

class VerletListMixedTabulatedLocal(InteractionLocal, interaction_VerletListMixedTabulated):
    def __init__(self, vl):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, interaction_VerletListMixedTabulated, vl)

    def setPotential(self, type1, type2, potential):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setPotential(self, type1, type2, potential)

    def getPotential(self, type1, type2):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            return self.cxxclass.getPotential(self, type1, type2)

if pmi.isController:
    class MixedTabulated(Potential):
        'The MixedTabulated potential.'
        pmiproxydefs = dict(
            cls='espressopp.interaction.MixedTabulatedLocal',
            pmiproperty=('mix_value', )
        )

    class VerletListMixedTabulated(Interaction):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls='espressopp.interaction.VerletListMixedTabulatedLocal',
            pmicall=['setPotential', 'getPotential']
        )