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

"""
******************************************
**espressopp.analysis.ChemicalConversion**
******************************************

The object that computes the number of particles of given type.

.. function:: espressopp.analysis.ChemicalConversion(system, particle_type, total_count)

            :param system: The system object
            :type system: espressopp.System
            :param particle_type: The particle type.
            :type particle_type: int
            :param total_count: The total number of particles of given type.
            :type total_count: int
"""

from espressopp.esutil import cxxinit
from espressopp import pmi

from espressopp.analysis.Observable import *  # NOQA
from _espressopp import analysis_ChemicalConversion


class ChemicalConversionLocal(ObservableLocal, analysis_ChemicalConversion):
    """The (local) compute of conversion."""
    def __init__(self, system, particle_type, total_count):
        if pmi.workerIsActive():
            cxxinit(self, analysis_ChemicalConversion, system, particle_type, total_count)

if pmi.isController:
    class ChemicalConversion(Observable):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls='espressopp.analysis.ChemicalConversionLocal',
            pmiproperty=['value']
        )
