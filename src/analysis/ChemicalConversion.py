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

Computes the ration of particles of given type with respect to the total number of particles.

.. math::

   x = n_t / n_{total}


where :math:`n_t` is a number of particles of type :math:`t` and :math:`n_{total}` is an input
parameter that in principle is a total number of particles of given type.

The usage is within the chemical reaction framework when during the reaction, the particles
change the types.

Two kinds of observables are implemented, basic one that counts the number of particles of given
type and second based on the sequence of types.

.. function:: espressopp.analysis.ChemicalConversion(system, particle_type, total_count)

            :param system: The system object
            :type system: espressopp.System
            :param particle_type: The particle type.
            :type particle_type: int
            :param total_count: The total number of particles of given type.
            :type total_count: int

.. function:: espressopp.analysis.ChemicalConversionTypeSequence(system, particle_group, type_seq, total_count)

            :param system: The system object
            :type system: espressopp.System
            :param particle_grouo: The particle type.
            :type particle_type: espressopp.ParticleGroup
            :param type_seq: The sequence of types to look for.
            :type type_seq: list
            :param total_count: The total number of molecules of given type sequence.
            :type total_count: int


Usage of ChemicalConversionTypeSequence
========================================

>>>> p = espressopp.ParticleGroup(system.storage)
>>>> [p.add(x) for x in range(1, 100)]   # creates a group of particles with ids 1...99
>>>> observable = espressopp.analysis.ChemicalConversionTypeSequence(system, p, [1, 2, 3], 100)

"""

from espressopp.esutil import cxxinit
from espressopp import pmi

from espressopp.analysis.Observable import *  # NOQA
from _espressopp import analysis_ChemicalConversion
from _espressopp import analysis_ChemicalConversionTypeSequence


class ChemicalConversionLocal(ObservableLocal, analysis_ChemicalConversion):
    """The (local) compute of conversion."""
    def __init__(self, system, particle_type, total_count):
        if pmi.workerIsActive():
            cxxinit(self, analysis_ChemicalConversion, system, particle_type, total_count)

class ChemicalConversionTypeSequenceLocal(ObservableLocal, analysis_ChemicalConversionTypeSequence):
    """The (local) compute of conversion."""
    def __init__(self, system, particle_group, type_sequence, total_count):
        if pmi.workerIsActive():
            cxxinit(self, analysis_ChemicalConversion, system, particle_group, total_count)

            self.cxxclass.set_sequence(type_sequence)

if pmi.isController:
    class ChemicalConversion(Observable):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls='espressopp.analysis.ChemicalConversionLocal',
        )

    class ChemicalConversionTypeSequence(Observable):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls='espressopp.analysis.ChemicalConversionTypeSequenceLocal',
            pmicall=['set_sequence']
        )
