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
******************************************
**espressopp.storage.DomainDecompositionFree**
******************************************


.. function:: espressopp.storage.DomainDecompositionFree(system, nodeGrid, cellGrid)

		:param system: 
		:param nodeGrid: 
		:param cellGrid: 
		:type system: 
		:type nodeGrid: 
		:type cellGrid: 

.. function:: espressopp.storage.DomainDecompositionFree.getCellGrid()

		:rtype: 

.. function:: espressopp.storage.DomainDecompositionFree.getNodeGrid()

		:rtype: 
"""
from espressopp import pmi
from espressopp.esutil import cxxinit
from _espressopp import storage_DomainDecompositionFree
from _espressopp import storage_DomainDecompositionFree
from espressopp import Int3D, toInt3DFromVector
from espressopp.tools import decomp
from espressopp import check
import mpi4py.MPI as MPI

from espressopp.storage.Storage import *

class DomainDecompositionFreeLocal(StorageLocal, storage_DomainDecompositionFree):

    def __init__(self, system, nodeGrid, cellGrid):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, storage_DomainDecompositionFree, system, nodeGrid, cellGrid)

    def getCellGrid(self):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            return self.cxxclass.getCellGrid(self)

    def getNodeGrid(self):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            return self.cxxclass.getNodeGrid(self)


if pmi.isController:
    class DomainDecompositionFree(Storage):
        pmiproxydefs = dict(
          cls = 'espressopp.storage.DomainDecompositionFreeLocal',
          pmicall = ['getCellGrid', 'getNodeGrid', 'cellAdjust']
        )
        def __init__(self, system, 
                     nodeGrid='auto', 
                     cellGrid='auto',
                     nocheck=False):
            # do sanity checks for the system first
            if nocheck:
              self.next_id = 0
              self.pmiinit(system, nodeGrid, cellGrid)                
            else:
              if check.System(system, 'bc'):
                if nodeGrid == 'auto':
                  nodeGrid = decomp.nodeGrid(system.comm.rank)
                else:
                  nodeGrid = toInt3DFromVector(nodeGrid)
                if cellGrid == 'auto':
                  cellGrid = Int3D(2,2,2)
                else:
                  cellGrid = toInt3DFromVector(cellGrid)
                # minimum image convention check:
                for k in range(3):
                  if nodeGrid[k]*cellGrid[k] == 1 :
                    print(("Warning! cellGrid[{}] has been "
                           "adjusted to 2 (was={})".format(k, cellGrid[k])))
                    cellGrid[k] = 2
                self.next_id = 0
                self.pmiinit(system, nodeGrid, cellGrid)
              else:
                print 'Error: could not create DomainDecomposition object'