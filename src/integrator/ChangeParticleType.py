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
********************************************
**espressopp.integrator.ChangeParticleType**
********************************************

This extension allows chaining particle type of selected number of particles during the MD simulation.
The change is done every *interval* steps for maximum *num_per_interval* particles.

.. function:: espressopp.integrator.ChangeParticleType(system, interval, num_per_interval, old_type_id, new_type_id)

		:param system: The system object.
		:param interval: Run every n-th MD steps.
		:type system: espressopp.System
		:type interval: int
		:param num_per_interval: Modify i particles.
		:type num_per_interval: int
		:param old_type_id: Particle type id to modify.
		:type old_type_id: int
		:param new_type_id: New particle type id.
		:type new_type_id: int

"""
from espressopp.esutil import cxxinit
from espressopp import pmi
from espressopp.integrator.Extension import *
from _espressopp import integrator_ChangeParticleType


class ChangeParticleTypeLocal(ExtensionLocal, integrator_ChangeParticleType):
    def __init__(self, system, interval, num_per_interval, old_type_id, new_type_id):
        if pmi.workerIsActive():
            cxxinit(self, integrator_ChangeParticleType, system, interval, num_per_interval, old_type_id, new_type_id)


if pmi.isController :
    class ChangeParticleType(Extension):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls =  'espressopp.integrator.ChangeParticleTypeLocal',
            pmicall = ['update_particles'])