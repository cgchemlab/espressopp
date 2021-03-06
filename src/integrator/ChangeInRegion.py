#  Copyright (C) 2016-2017
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
****************************************
**espressopp.integrator.ChangeInRegion**
****************************************

This extension allows chainging properties of particles whenever they enters the defined geometrical region.

.. function:: espressopp.integrator.ChangeInRegion(system, particleGroup, fixMask)

		:param system: 
		:param particleGroup: 
		:type system:
		:type particleGroup:

.. function:: espressopp.integrator.ChangeInRegion.set_particle_properties(type_id, particle_properties)

      Defines the change of the particle properties, the particles of `type_id` in region will be updated only.

      :param type_id: The particle type id.
      :type type_id:

.. function:: espressopp.integrator.ChangeInRegion.set_flags(type_id, reset_velocity, reset_force, remove_particle)

      Sets if velocity of force of particle of given type should have reset the force and velocity.

      :param type_id: The particle type id.
      :type type_id: int
      :param reset_velocity: If set to true then the velocity will be reset to 0.0.
      :type reset_velocity: bool
      :param reset_force: If set to true the the force will be reset to 0.0.
      :type reset_force: bool
      :param remove_particle: If set to true then the particle will be removed from the simulation
      :type remove_particle: bool

"""
from espressopp.esutil import cxxinit
from espressopp import pmi
from espressopp.integrator.Extension import *
from _espressopp import integrator_ChangeInRegion

class ChangeInRegionLocal(ExtensionLocal, integrator_ChangeInRegion):
    def __init__(self, system, particleGroup, prob=None, p_num=None, p_num_percentage=None):
        if pmi.workerIsActive():
            if prob:
                cxxinit(self, integrator_ChangeInRegion, system, particleGroup, prob)
            elif p_num:
                cxxinit(self, integrator_ChangeInRegion, system, particleGroup, p_num)
            elif p_num_percentage:
                cxxinit(self, integrator_ChangeInRegion, system, particleGroup, 0, p_num_percentage)
            else:
                raise RuntimeError('prob, p_num, p_percentage have to be specified')

    def set_flags(self, type_id, reset_velocity, reset_force, remove_particle=False):
        if pmi.workerIsActive():
            self.cxxclass.set_flags(self, type_id, reset_velocity, reset_force, remove_particle)

    def set_particle_properties(self, type_id, particle_properties):
        if pmi.workerIsActive():
            self.cxxclass.set_particle_properties(self, type_id, particle_properties)

if pmi.isController :
    class ChangeInRegion(Extension):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls =  'espressopp.integrator.ChangeInRegionLocal',
            pmicall = ['set_particle_properties', 'update_particles', 'set_flags'],
            pmiproperty = ('stats_filename',)
            )
