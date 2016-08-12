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
**************************************
**FreeOrthorhombicBC** - Object
**************************************

Like all boundary condition objects, this class implements
all the methods of the base class **BC** , which are described in detail
in the documentation of the abstract class **BC**.

The FreeOrthorhombicBC class is responsible for the FreeOrthorhombic boundary condition.
Currently only periodic boundary conditions are supported.

Example: 

>>> boxsize = (Lx, Ly, Lz)
>>> bc = espressopp.bc.FreeOrthorhombicBC(rng, boxsize)


.. function:: espressopp.bc.FreeOrthorhombicBC(rng, boxL, periodicity)

		:param rng: 
		:param boxL: (default: 1.0)
		:type rng: 
		:type boxL: real

.. function:: espressopp.bc.FreeOrthorhombicBC.setBoxL(boxL)

		:param boxL: 
		:type boxL: 
"""

from espressopp.esutil import cxxinit
from espressopp import pmi
from espressopp import toReal3D
from espressopp import toInt3D

from espressopp.bc.BC import *
from _espressopp import bc_FreeOrthorhombicBC

class FreeOrthorhombicBCLocal(BCLocal, bc_FreeOrthorhombicBC):
    def __init__(self, rng, boxL=1.0, periodicity=[True, True, True]):
        if pmi.workerIsActive():
            cxxinit(self, bc_FreeOrthorhombicBC, rng, toReal3D(boxL), toInt3D(periodicity))

    # override length property
    def setBoxL(self, boxL):
        if pmi.workerIsActive():
            self.cxxclass.boxL.fset(self, toReal3D(boxL))

    boxL = property(bc_FreeOrthorhombicBC.boxL.fget, setBoxL)

if pmi.isController :
    class FreeOrthorhombicBC(BC):
        pmiproxydefs = dict(
            cls = 'espressopp.bc.FreeOrthorhombicBCLocal',
            pmiproperty = ['boxL']
            )
