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
*************************************
espressopp.analysis.AngleDistribution
*************************************

Calculates angle distribution on the fly during the simulation, using the data stored
in the fixed triple list.


.. function:: espressopp.analysis.AngleDistribution(system)

		:param system:
		:type system:

.. function:: espressopp.analysis.AngleDistribution.compute(n)

        :param n: number of bins in the histogram
        :type n: int


"""
from espressopp.esutil import cxxinit
from espressopp import pmi

from espressopp.analysis.Observable import *
from _espressopp import analysis_AngleDistribution


class AngleDistributionLocal(ObservableLocal, analysis_AngleDistribution):
    def __init__(self, system):
        if pmi.workerIsActive():
            cxxinit(self, analysis_AngleDistribution, system)

    def compute(self, n):
        if pmi.workerIsActive():
            return self.cxxclass.compute(self, n)

    def load_from_topology_manager(self, tm):
        if pmi.workerIsActive():
            self.cxxclass.load_from_topology_manager(self, tm)

if pmi.isController:
    class AngleDistribution(Observable):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            pmicall=["compute", "register_triplet", "load_from_topology_manager"],
            cls='espressopp.analysis.AngleDistributionLocal'
        )
