#
#  Copyright (C) 2015-2017
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
*********************************
**espressopp.integrator.TopologyManager**
*********************************


.. function:: espressopp.integrator.TopologyManager(system)

		:param system: 
		:param verletlist: 
		:type system: 
		:type verletlist: 
"""
from espressopp.esutil import cxxinit
from espressopp import pmi

from _espressopp import integrator_TopologyManager
from _espressopp import integrator_TopologyParticleProperties


class TopologyParticlePropertiesLocal(integrator_TopologyParticleProperties):
    def __init__(self, type=None, mass=None, q=None, lambda_adr=None, incr_state=None, state=None):
        if pmi.workerIsActive():
            cxxinit(self, integrator_TopologyParticleProperties)
            if incr_state is not None and state is not None:
                raise RuntimeError('Ambiguity, cannot set incr_state and state at the same time')
            if type is not None:
                self.type_id = int(type)
            if mass is not None:
                self.mass = mass
            if q is not None:
                self.q = q
            if lambda_adr is not None:
                self.lambda_adr = lambda_adr
            if incr_state is not None:
                self.incr_state = incr_state
            if state is not None:
                self.state = state

    def set_min_max_state(self, min_state, max_state):
        if pmi.workerIsActive():
            self.cxxclass.set_min_max_state(self, min_state, max_state)


class TopologyManagerLocal(integrator_TopologyManager):

    def __init__(self, system):
        if pmi.workerIsActive():
            cxxinit(self, integrator_TopologyManager, system)

    def observe_tuple(self, fpl):
        if pmi.workerIsActive():
            self.cxxclass.observe_tuple(self, fpl)

    def register_tuple(self, fpl, type1, type2=None):
        if pmi.workerIsActive():
            if type2 is None:
                type2 = type1
            self.cxxclass.register_tuple(self, fpl, type1, type2)

    def register_14tuple(self, fpl, type1, type2=None):
        if pmi.workerIsActive():
            if type2 is None:
                type2 = type1
            self.cxxclass.register_14tuple(self, fpl, type1, type2)

    def register_triplet(self, ftl, type1, type2=None, type3=None):
        if pmi.workerIsActive():
            if type2 is None and type3 is None:
                type2 = type3 = type1
            self.cxxclass.register_triple(self, ftl, type1, type2, type3)

    def register_quadruplet(self, fql, type1, type2=None, type3=None, type4=None):
        if pmi.workerIsActive():
            if None in [type2, type3, type4]:
                type2 = type3 = type4 = type1
            self.cxxclass.register_quadruple(self, fql, type1, type2, type3, type4)

    def initialize_topology(self):
        if pmi.workerIsActive():
            self.cxxclass.initialize(self)

    def is_residue_connected(self, rid1, rid2):
        if pmi.workerIsActive():
            return self.cxxclass.is_residue_connected(self, rid1, rid2)

    def is_particle_connected(self, pid1, pid2):
        if pmi.workerIsActive():
            return self.cxxclass.is_particle_connected(self, pid1, pid2)

    def get_fixed_pair_list(self, type1, type2):
        if pmi.workerIsActive():
            return self.cxxclass.get_fixed_pair_list(self, type1, type2)


if pmi.isController :
    class TopologyParticleProperties(object):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls='espressopp.integrator.TopologyParticlePropertiesLocal',
            pmiproperty=['type_id', 'mass', 'q', 'state', 'lambda_adr', 'incr_state'],
            pmicall=('set_min_max_state',)
        )

        def __str__(self):
            obj = self.pmiobject
            return 'ParticleProperties(type={}, mass={}, q={}, lambda={}, incr_state={})'.format(
                obj.type_id, obj.mass, obj.q, obj.lambda_adr, obj.incr_state)

    class TopologyManager(object):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls =  'espressopp.integrator.TopologyManagerLocal',
            pmicall = ['observe_tuple', 'register_tuple', 'register_14tuple', 'register_triplet',
                       'register_quadruplet', 'initialize_topology', 'exchange_data',
                       'is_residue_connected', 'is_particle_connected',
                       'save_topology', 'save_res_topology', 'save_residues',
                       'has_neighbour_particle_property',
                       'get_fixed_pair_list', 'get_fixed_triple_list'
                      ],
            pmiinvoke = [
                'print_topology',
                'print_res_topology',
                'print_residues',
                'get_neighbour_lists',
                'get_timers',
                'get_molecule_ids',
                'get_molecule',
                'get_residue_id',
                'get_molecule_id']
            )
