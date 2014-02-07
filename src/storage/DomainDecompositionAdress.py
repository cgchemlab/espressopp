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
**************************************
**DomainDecompositionAdress** - Object
**************************************

The DomainDecompositionAdress is the Domain Decomposition for AdResS and H-
AdResS simulations. It makes sure that tuples (i.e. a coarse-grained particle
and its corresponding atomistic particles) are always stored together on one CPU.
When setting DomainDecompositionAdress you have to provide the system as well as
the nodegrid and the cellgrid.

Example - setting DomainDecompositionAdress:

>>> system.storage = espresso.storage.DomainDecompositionAdress(system, nodeGrid, cellGrid)

"""

from espresso import pmi
from espresso.esutil import cxxinit
from _espresso import storage_DomainDecompositionAdress
from espresso import Int3D, toInt3DFromVector
import mpi4py.MPI as MPI

from espresso.storage.Storage import *

class DomainDecompositionAdressLocal(StorageLocal, 
                               storage_DomainDecompositionAdress):
    'The (local) DomainDecomposition.'
    def __init__(self, system, nodeGrid, cellGrid):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, storage_DomainDecompositionAdress, system, nodeGrid, cellGrid)

if pmi.isController:
    class DomainDecompositionAdress(Storage):
        pmiproxydefs = dict(
            cls = 'espresso.storage.DomainDecompositionAdressLocal',
            pmicall = ['getCellGrid', 'cellAdjust']
            )
        def __init__(self, system, 
                     nodeGrid='auto', 
                     cellGrid='auto'):
            if nodeGrid == 'auto':
                nodeGrid = Int3D(system.comm.rank, 1, 1)
            else:
                nodeGrid = toInt3DFromVector(nodeGrid)

            if cellGrid == 'auto':
                # TODO: Implement
                raise 'Automatic cell size calculation not yet implemented'
            else:
                cellGrid = toInt3DFromVector(cellGrid)

            self.next_id = 0
            self.pmiinit(system, nodeGrid, cellGrid)
