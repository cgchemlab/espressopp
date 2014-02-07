#  Copyright (C) 2012,2013
#      Max Planck Institute for Polymer Research
#  Copyright (C) 2008,2009,2010,2011
#      Max-Planck-Institute for Polymer Research & Fraunhofer SCAI
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
****************************************
**espresso.storage.DomainDecomposition**
****************************************

"""
from espresso import pmi
from espresso.esutil import cxxinit
from _espresso import storage_DomainDecomposition
from espresso import Int3D, toInt3DFromVector
from espresso.tools import decomp
from espresso import check
import mpi4py.MPI as MPI

from espresso.storage.Storage import *

class DomainDecompositionLocal(StorageLocal, storage_DomainDecomposition):
    'The (local) DomainDecomposition.'
    def __init__(self, system, nodeGrid, cellGrid):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, storage_DomainDecomposition, system, nodeGrid, cellGrid)
    
    def getCellGrid(self):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            return self.cxxclass.getCellGrid(self)

    def getNodeGrid(self):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            return self.cxxclass.getNodeGrid(self)
          
if pmi.isController:
    class DomainDecomposition(Storage):
        pmiproxydefs = dict(
          cls = 'espresso.storage.DomainDecompositionLocal',  
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
                    cellGrid[k] = 2
                    print 'cellGrid[%i] has been adjusted to 2'                  
                self.next_id = 0
                self.pmiinit(system, nodeGrid, cellGrid)
              else:
                print 'Error: could not create DomainDecomposition object'
