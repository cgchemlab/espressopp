#  Copyright (C) 2015
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
***************************************
**espressopp.analysis.NFixedPairListEntries**
***************************************

The object that computes the number of entries in FixedPairList.

.. function:: espressopp.analysis.NFixedPairListEntries(system, fixed_pair_lit)

            :param system: The system object
            :type system: espressopp.System
            :param fixed_pair_list: The observed fpl.
            :type interaction: espressopp.FixedPairList
"""

from espressopp.esutil import cxxinit
from espressopp import pmi

from espressopp.analysis.Observable import *  # NOQA
from _espressopp import analysis_NFixedPairListEntries
from _espressopp import analysis_NFixedTripleListEntries
from _espressopp import analysis_NFixedQuadrupleListEntries


class NFixedPairListEntriesLocal(ObservableLocal, analysis_NFixedPairListEntries):
    """The (local) compute of potential energy."""
    def __init__(self, system, fl):
        if pmi.workerIsActive():
            cxxinit(self, analysis_NFixedPairListEntries, system, fl)

class NFixedTripleListEntriesLocal(ObservableLocal, analysis_NFixedTripleListEntries):
    """The (local) compute of potential energy."""
    def __init__(self, system, fl):
        if pmi.workerIsActive():
            cxxinit(self, analysis_NFixedTripleListEntries, system, fl)

class NFixedQuadrupleListEntriesLocal(ObservableLocal, analysis_NFixedQuadrupleListEntries):
    """The (local) compute of potential energy."""
    def __init__(self, system, fl):
        if pmi.workerIsActive():
            cxxinit(self, analysis_NFixedQuadrupleListEntries, system, fl)

if pmi.isController:
    class NFixedPairListEntries(Observable):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls='espressopp.analysis.NFixedPairListEntriesLocal',
            pmiproperty=['value']
        )

    class NFixedTripleListEntries(Observable):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls='espressopp.analysis.NFixedTripleListEntriesLocal',
            pmiproperty=['value']
        )

    class NFixedQuadrupleListEntries(Observable):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls='espressopp.analysis.NFixedQuadrupleListEntriesLocal',
            pmiproperty=['value']
        )
