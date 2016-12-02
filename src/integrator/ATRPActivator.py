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


"""
***************************************
**espressopp.integrator.ATRPActivator**
***************************************

This extension allows activating or deactivating with given probability particles of given type in given chemical
state.



.. function:: espressopp.integrator.ATRPActivator(system, interval, num_per_interval, ratio_activator,
                                                  ratio_deactivator, delta_catalyst, k_activate,
                                                  k_deactivate, old_type_id, new_type_id)

		:param system: The system object.
		:param interval: Run every n-th MD steps.
		:type system: espressopp.System
		:type interval: int
		:param num_per_interval: Modify i-th particles per interval.
		:type num_per_interval: int
		:param ratio_activator: ratio of the CuI to total number of Cu.
		:type ratio_activator:float
		:param ratio_deactivator: ratio of the CuII to total number of Cu.
		:type ratio_deactivator:float
		:param delta_catalyst: the increasement of ratio_activator/deactivator.
		:type delta_catalyst: float
		:param k_activate: activate constant.
		:type k_activate: float
		:param k_deactivate: deactivate constant.
		:type k_deactivate: float

.. function:: espressopp.integrator.ATRPActivator.add_reactive_center(type_id, state, is_activator,
                                                                      new_property, delta_state)
        Add defintion of reactive center and change of properties.

        :param type_id: The particle type.
        :type type_id: int
        :param state: The state.
        :type state: int
        :param is_activator:
        :type is_activator: bool
        :param new_property: The new property of particle.
        :type new_property: espressopp.integrator.TopologyParticleProperties
        :param delta_state: The change of chemical state.
        :type delta_state: int

"""
from espressopp.esutil import cxxinit
from espressopp import pmi
from espressopp.integrator.Extension import *
from _espressopp import integrator_ATRPActivator


class ATRPActivatorLocal(ExtensionLocal, integrator_ATRPActivator):
    def __init__(self, system, interval, num_per_interval, ratio_activator,
                 ratio_deactivator, delta_catalyst, k_activate, k_deactivate):
        if pmi.workerIsActive():
            cxxinit(self, integrator_ATRPActivator, system, interval, num_per_interval, ratio_activator,
                    ratio_deactivator, delta_catalyst, k_activate, k_deactivate)

    def add_reactive_center(self, type_id, state, is_activator, new_property, delta_state):
        """Defines reactive center"""
        if pmi.workerIsActive():
            self.cxxclass.add_reactive_center(
                self, type_id, state, is_activator, new_property, delta_state)

if pmi.isController :
    class ATRPActivator(Extension):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls =  'espressopp.integrator.ATRPActivatorLocal',
            pmiproperty = ['stats_filename'],
            pmicall = ['update_particles', 'add_reactive_center', 'save_statistics'])