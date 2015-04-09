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
************************************
**FixedPairListType** - Object
************************************

FixedPairListType is a special case of standard FixedPairList. It keeps
only bonds for particles of defined types. In this way, whenever particle
change the type during simulation, the bond will be removed if it is not valid.

Example - creating the FixedPairListType and adding bonds:

>>> fpl = espressopp.FixedPairListType(system.storage, 1, 2)
>>> fpl.addBonds(bonds)
"""

from espressopp import pmi
import _espressopp
import espressopp
from espressopp.esutil import cxxinit


class FixedPairListTypeLocal(_espressopp.FixedPairListType):
    'The (local) fixed pair list.'

    def __init__(self, storage, type1, type2):
        'Local construction of a fixed pair list'
        if pmi.workerIsActive():
            cxxinit(self, _espressopp.FixedPairListType, storage, type1, type2)

    def add(self, pid1, pid2):
        'add pair to fixed pair list'
        if pmi.workerIsActive():
            return self.cxxclass.add(self, pid1, pid2)

    def getBonds(self):
        'return the bonds of the GlobalPairList'
        if pmi.workerIsActive():
          bonds=self.cxxclass.getBonds(self)
          return bonds

    def addBonds(self, bondlist):
        """
        Each processor takes the broadcasted bondlist and
        adds those pairs whose first particle is owned by
        this processor.
        """
        if pmi.workerIsActive():
            for bond in bondlist:
                pid1, pid2 = bond
                self.cxxclass.add(self, pid1, pid2)


if pmi.isController:
    class FixedPairListType(object):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls = 'espressopp.FixedPairListTypeLocal',
            localcall = [ "add" ],
            pmicall = [ "addBonds" ],
            pmiinvoke = ['getBonds']
            )
