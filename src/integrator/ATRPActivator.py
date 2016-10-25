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



.. function:: espressopp.integrator.ATRPActivator(system, interval, num_per_interval, old_type_id, new_type_id)

		:param system: The system object.
		:param interval: Run every n-th MD steps.
		:type system: espressopp.System
		:type interval: int
		:param num_per_interval: Modify i-th particles per interval.
		:type num_per_interval: int

.. function:: espressopp.integrator.ATRPActivator.add_reactive_center(type_id, min_state, max_state,
                                                                      new_property, delta_state, prob)
        Add defintion of reactive center and change of properties.

        :param type_id: The particle type.
        :type type_id: int
        :param min_state: The min state.
        :type min_state: int
        :param max_state: The max state.
        :type max_state: int
        :param new_property: The new property of particle.
        :type new_property: espressopp.ParticleProperties
        :param delta_state: The change of chemical state.
        :type delta_state: int
        :param prob: The probability constant.
        :type prob: float

"""
from espressopp.esutil import cxxinit
from espressopp import pmi
from espressopp.integrator.Extension import *
from _espressopp import integrator_ATRPActivator


class ATRPActivatorLocal(ExtensionLocal, integrator_ATRPActivator):
    def __init__(self, system, interval, num_per_interval):
        if pmi.workerIsActive():
            cxxinit(self, integrator_ATRPActivator, system, interval, num_per_interval)

    def add_reactive_center(self, type_id, min_state, max_state, new_property, delta_state, prob):
        """Defines reactive center"""
        if pmi.workerIsActive():
            self.cxxclass.add_reactive_center(
                self, type_id, min_state, max_state, new_property, delta_state, prob)

if pmi.isController :
    class ATRPActivator(Extension):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls =  'espressopp.integrator.ATRPActivatorLocal',
            pmicall = ['update_particles', 'add_reactive_center'])