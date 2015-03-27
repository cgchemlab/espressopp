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
********************************************
**espressopp.interaction.LennardJonesGromacs**
********************************************

"""
from espressopp import pmi, infinity
from espressopp.esutil import *

from espressopp.interaction.Potential import *
from espressopp.interaction.Interaction import *
from _espressopp import interaction_LennardJonesGromacs, \
                      interaction_VerletListLennardJonesGromacs, \
                      interaction_CellListLennardJonesGromacs, \
                      interaction_FixedPairListLennardJonesGromacs

class LennardJonesGromacsLocal(PotentialLocal, interaction_LennardJonesGromacs):
    'The (local) LennardJonesGromacs potential.'
    def __init__(self, epsilon=1.0, sigma=1.0, r1=0.0,
                 cutoff=infinity, shift="auto"):
        """Initialize the local LennardJonesGromacs object."""
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            if shift =="auto":
                cxxinit(self, interaction_LennardJonesGromacs,
                        epsilon, sigma, r1, cutoff)
            else:
                cxxinit(self, interaction_LennardJonesGromacs,
                        epsilon, sigma, r1, cutoff, shift)

class VerletListLennardJonesGromacsLocal(InteractionLocal, interaction_VerletListLennardJonesGromacs):
    'The (local) LennardJonesGromacs interaction using Verlet lists.'
    def __init__(self, vl):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, interaction_VerletListLennardJonesGromacs, vl)

    def setPotential(self, type1, type2, potential):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setPotential(self, type1, type2, potential)

    def getPotential(self, type1, type2):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            return self.cxxclass.getPotential(self, type1, type2)

class CellListLennardJonesGromacsLocal(InteractionLocal, interaction_CellListLennardJonesGromacs):
    'The (local) LennardJonesGromacs interaction using cell lists.'
    def __init__(self, stor):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, interaction_CellListLennardJonesGromacs, stor)
        
    def setPotential(self, type1, type2, potential):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setPotential(self, type1, type2, potential)

class FixedPairListLennardJonesGromacsLocal(InteractionLocal, interaction_FixedPairListLennardJonesGromacs):
    'The (local) LennardJonesGromacs interaction using FixedPair lists.'
    def __init__(self, system, vl, potential):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, interaction_FixedPairListLennardJonesGromacs, system, vl, potential)
        
    def setPotential(self, potential):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setPotential(self, potential)

if pmi.isController:
    class LennardJonesGromacs(Potential):
        'The LennardJonesGromacs potential.'
        pmiproxydefs = dict(
            cls = 'espressopp.interaction.LennardJonesGromacsLocal',
            pmiproperty = ['epsilon', 'sigma', 'r1']
            )

    class VerletListLennardJonesGromacs(Interaction):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls =  'espressopp.interaction.VerletListLennardJonesGromacsLocal',
            pmicall = ['setPotential','getPotential']
            )
    class CellListLennardJonesGromacs(Interaction):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls =  'espressopp.interaction.CellListLennardJonesGromacsLocal',
            pmicall = ['setPotential']
            )
    class FixedPairListLennardJonesGromacs(Interaction):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls =  'espressopp.interaction.FixedPairListLennardJonesGromacsLocal',
            pmicall = ['setPotential']
            )
