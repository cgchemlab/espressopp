#  Copyright (C) 2014-2016
#    Jakub Krajniak (jkrajniak@gmail.com)
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
##########################################
**Integrator extension**
##########################################

.. function:: espressopp.integrator.ChemicalReaction(system, vl, domdec, interval)

    :param system: The system object.
    :param vl: The verlet list object.
    :param domdec: The domain decomposition object.
    :param interval: The time between reactions are invoked. (:math:`\Phi`)
    :type system: espressopp.System
    :type vl: espressopp.VerletList
    :type domdec: espressopp.storage.DomainDecomposition
    :type interval: integer

.. function:: espressopp.integrator.ChemicaReaction.add_reaction(reaction):

        Adds chemical reaction.

    :param reaction: The Reaction object.
    :type reaction: espressopp.integrator.Reaction

==================================
**espressopp.integrator.Reaction**
==================================

.. function:: espressopp.integrator.Reaction(type_1, type_2, delta_1, delta_2, min_state_1, max_state_1,
                                             min_state_2, max_state_2, cutoff, rate, fpl)

    Definition of chemical reaction that follows the equation:

    .. math::   A_{T1}^{i} + B_{T2}^{j} \rightarrow A^{i+\delta i}:B^{j+\delta j}

    with the default distance criterion:

    .. math:: \|\mathbf{r}_A - \mathbf{r}_B\| < r_c

    and

    .. math:: i_{min} \le i < i_{max} \land j_{min} \le j < j_{max}

    The reaction rate is given with `rate` parameter. It is possible to define rate parameter that
    will depend on the chemical state of particle A or B. The rate can be defined separately for both
    particle types by passing a tuple of two lists

        ([i_min, i_min+1, i_min+2, ... , i_max-1], [j_min, j_min+1, j_min+2, ... , j_max-1])

    If rate is defined as a single positive real number then it will be assigned to all states, both for A and
    B particles.

    :param type_1: The type of particle A.
    :type type_1: int
    :param type_2: The type of particle B.
    :type type_2: int
    :param delta_1: The change on state of particle A.
    :type delta_1: int
    :param delta_2: The change on state of particle A.
    :type delta_2: int
    :param min_state_1: The minimum of state A.
    :type min_state_1: int
    :param max_state_1: The maximum of state A.
    :type max_state_1: int
    :param min_state_2: The minimum of state B.
    :type max_state_2: The maximum of state B.
    :param cutoff: The cutoff distance.
    :type cutoff: float
    :param rate: The reaction rate.
    :type rate: float
    :param fpl: The fixed tuple list that will store new bonds.
    :type fpl: espressopp.FixedPairList


.. function:: espressopp.integrator.Reaction.add_postprocess(post_process, type)

    Adds post-process object that will be invoked on every new pair formed by
    this reaction. The type defines on which particle in the pair A, B, the post-process
    object will be invoked:

      - `both`: On both, A and B particles.
      - `type_1`: On A particle (defined type `type_1` in the reaction)
      - `type_2`: On B particle (defined type `type_2` in the reaction)

    :param post_process: The post process object.
    :type post_process: espressopp.integrator.ChemicalReactionPostProcess
    :param type: Defined on which particle post-process will be run.
    :type type: string

.. function:: espressopp.integrator.Reaction.set_reaction_cutoff(rc)

    Define how the distance condition will be handled.

    :param rc: The ReactionCutoff object.
    :type rc: espressopp.integrator.ReactionCutoffStatic or espressopp.integrator.ReactionCutoffRandom


========================================
**Reaction cut-off**
========================================

This defines how the distance criterion will be handled. There are two ways, one is static
and second is based on probability distrbution function (see image: :ref:`fig_distance_condition`)

There are two classes:

 - **espressopp.integrator.ReactionCutoffStatic**
 - **espressopp.integrator.ReactionCutoffRandom**

++++++++++++++++++++++++++++++++++++++++++++++
**espressopp.integrator.ReactionCutoffStatic**
++++++++++++++++++++++++++++++++++++++++++++++

.. function:: espressopp.integrator.ReactionCutoffStatic(max_cutoff, min_cutoff=0.0)

   :param max_cutoff: The :math:`r_2` distance.
   :type max_cutoff: float
   :param min_cutoff: The :math:`r_1` distance.
   :type min_cutoff: float


++++++++++++++++++++++++++++++++++++++++++++++
**espressopp.integrator.ReactionCutoffRandom**
++++++++++++++++++++++++++++++++++++++++++++++

Whenever :ref:`espressopp.integrator.Reaction` checks the distance condition between two particles,
a random number `W` is selected from normal distribution with :math:`\mu = 0.0` and :math:`\sigma = eq\_width`
where `eq_width` is a width of the distribution.

The condition is as follows:

  .. math:: \|\mathbf{r}_A - \mathbf{r}_B\| < |W| + b_0

.. function:: espressopp.integrator.ReactionCutoffRandom(eq_distance, eq_width, seed=12345)

   :param eq_distance: The :math:`b_0` value.
   :type eq_distance: float
   :param eq_width: The :math:`eq\_width` value.
   :type eq_width: float
"""


from espressopp.esutil import cxxinit
from espressopp import pmi

from espressopp.integrator.Extension import *  # NOQA
from _espressopp import integrator_ChemicalReaction
from _espressopp import integrator_Reaction
from _espressopp import integrator_RestrictReaction
from _espressopp import integrator_DissociationReaction

from _espressopp import integrator_ChemicalReactionPostProcess
from _espressopp import integrator_PostProcessChangeProperty
from _espressopp import integrator_PostProcessRemoveBond
from _espressopp import integrator_PostProcessChangeNeighboursProperty
from _espressopp import integrator_PostProcessChangePropertyOnState
from _espressopp import integrator_PostProcessRemoveNeighbourBond

from _espressopp import integrator_ReactionCutoffStatic
from _espressopp import integrator_ReactionCutoffRandom
from _espressopp import integrator_ReactionCutoff


class ChemicalReactionLocal(ExtensionLocal, integrator_ChemicalReaction):
    """Chemical Reaction integrator extension."""
    def __init__(self, system, vl, domdec, topology_manager, interval=None):
        """Chemical reaction extension.

        Args:
          system: The espressopp.System object.
          vl: The verlet list.
          domdec: The domain decomposition object.
          topology_manager: The TopologyManager object.
          interval: The timestep where the extension will be run.
        """
        if pmi.workerIsActive():
            cxxinit(self, integrator_ChemicalReaction, system, vl, domdec, topology_manager)
            if interval is not None:
                self.interval = interval

    def add_reaction(self, reaction):
        """Adds the reaction to the list."""
        if pmi.workerIsActive():
            self.cxxclass.add_reaction(self, reaction)


class PostProcessChangePropertyLocal(integrator_PostProcessChangeProperty,
                                     integrator_ChemicalReactionPostProcess):
    """Post process of reaction that changes particle property."""
    def __init__(self, type_id=None, prop=None):
        if pmi.workerIsActive():
            cxxinit(self, integrator_PostProcessChangeProperty)
            if type_id is not None and prop is not None:
                self.add_change_property(type_id, prop)

    def add_change_property(self, type_id, prop):
        if pmi.workerIsActive():
            self.cxxclass.add_change_property(self, type_id, prop)

    def remove_change_property(self, type_id):
        if pmi.workerIsActive():
            self.cxxclass.remove_change_property(self, type_id)


class PostProcessChangeNeighboursPropertyLocal(integrator_PostProcessChangeNeighboursProperty,
                                               integrator_ChemicalReactionPostProcess):
    """Post process of reaction that changes particle property."""
    def __init__(self, topology_manager):
        if pmi.workerIsActive():
            cxxinit(self, integrator_PostProcessChangeNeighboursProperty, topology_manager)

    def add_change_property(self, type_id, prop, nb_level):
        if pmi.workerIsActive():
            self.cxxclass.add_change_property(self, type_id, prop, nb_level)


class PostProcessChangePropertyOnStateLocal(integrator_PostProcessChangePropertyOnState,
                                            integrator_ChemicalReactionPostProcess):
    """Post process of reaction that changes particle property."""
    def __init__(self):
        if pmi.workerIsActive():
            cxxinit(self, integrator_PostProcessChangePropertyOnState)

    def add_change_property(self, type_id, prop, nb_level):
        if pmi.workerIsActive():
            self.cxxclass.add_change_property(self, type_id, prop, nb_level)


class PostProcessRemoveNeighbourBondLocal(integrator_PostProcessRemoveNeighbourBond,
                                          integrator_ChemicalReactionPostProcess):
    def __init__(self, topology_manager):
        if pmi.workerIsActive():
            cxxinit(self, integrator_PostProcessRemoveNeighbourBond, topology_manager)

    def add_bond_to_remove(self, type_id, nb_level, type_pid1, type_pid2):
        if pmi.workerIsActive():
            self.cxxclass.add_bond_to_remove(self, type_id, nb_level, type_pid1, type_pid2)


class PostProcessRemoveBondLocal(integrator_PostProcessRemoveBond, integrator_ChemicalReactionPostProcess):
    """Post process of reaction, remove some bonds."""
    def __init__(self, fpl, number_of_bonds):
        if pmi.workerIsActive():
            cxxinit(self, integrator_PostProcessRemoveBond, fpl, number_of_bonds)


class ReactionCutoffStaticLocal(integrator_ReactionCutoffStatic, integrator_ReactionCutoff):
    """Reaction cutoff object."""
    def __init__(self, max_cutoff, min_cutoff=0.0):
        if pmi.workerIsActive():
            cxxinit(self, integrator_ReactionCutoffStatic, min_cutoff, max_cutoff)


class ReactionCutoffRandomLocal(integrator_ReactionCutoffRandom, integrator_ReactionCutoff):
    """Reaction cutoff from normal distribution."""
    def __init__(self, eq_distance, eq_width, seed=12345):
        if pmi.workerIsActive():
            cxxinit(self, integrator_ReactionCutoffRandom, eq_distance, eq_width, seed)


class ReactionLocal(integrator_Reaction):
    """Synthesis reaction."""
    def __init__(self, type_1, type_2, delta_1, delta_2, min_state_1, max_state_1,
                 min_state_2, max_state_2, cutoff, rate, fpl):
        if pmi.workerIsActive():
            cxxinit(
                self,
                integrator_Reaction,
                type_1,
                type_2,
                delta_1,
                delta_2,
                min_state_1,
                max_state_1,
                min_state_2,
                max_state_2,
                fpl,
                rate
            )
            self.cxxclass.set_reaction_cutoff(self, ReactionCutoffStaticLocal(cutoff))

    def add_postprocess(self, post_process, reactant_switch=0):
        """Add new post process to the reaction.
            Args:
                post_process: The post process object.
                reactant_switch: Which reactant (A, B or both) should be affected.
                    0 - both, 1 - reactant A, 2 - reactant B.
        """
        if pmi.workerIsActive():
            name_switch = {'both': 0, 'type_1': 1, 'type_2': 2}
            self.cxxclass.add_postprocess(self, post_process, name_switch.get(reactant_switch, reactant_switch))

    def set_reaction_cutoff(self, reaction_cutoff):
        """Set cutoff object."""
        if pmi.workerIsActive():
            self.cxxclass.set_reaction_cutoff(self, reaction_cutoff)


class RestrictReactionLocal(integrator_Reaction):
    """Synthesis reaction."""
    def __init__(self, type_1, type_2, delta_1, delta_2, min_state_1, max_state_1,
                 min_state_2, max_state_2, cutoff, rate, fpl):
        if pmi.workerIsActive():
            cxxinit(
                self,
                integrator_RestrictReaction,
                type_1,
                type_2,
                delta_1,
                delta_2,
                min_state_1,
                max_state_1,
                min_state_2,
                max_state_2,
                fpl,
                rate
            )
            self.cxxclass.set_reaction_cutoff(self, ReactionCutoffStaticLocal(cutoff))

    def add_postprocess(self, post_process, reactant_switch=0):
        """Add new post process to the reaction.
            Args:
                post_process: The post process object.
                reactant_switch: Which reactant (A, B or both) should be affected.
                    0 - both, 1 - reactant A, 2 - reactant B.
        """
        if pmi.workerIsActive():
            name_switch = {'both': 0, 'type_1': 1, 'type_2': 2}
            self.cxxclass.add_postprocess(self, post_process, name_switch.get(reactant_switch, reactant_switch))

    def set_reaction_cutoff(self, reaction_cutoff):
        """Set cutoff object."""
        if pmi.workerIsActive():
            self.cxxclass.set_reaction_cutoff(self, reaction_cutoff)

    def define_connection(self, pid1, pid2):
        if pmi.workerIsActive():
            self.cxxclass.define_connection(self, pid1, pid2)

    def define_connections(self, connection_list):
        if pmi.workerIsActive():
            for p1, p2 in connection_list:
                self.cxxclass.define_connection(self, p1, p2)


class DissociationReactionLocal(integrator_DissociationReaction):
    """DissociationReaction reaction."""
    def __init__(self, type_1, type_2, delta_1, delta_2, min_state_1, max_state_1,
                 min_state_2, max_state_2, cutoff, rate, fpl):
        if pmi.workerIsActive():
            cxxinit(
                self,
                integrator_DissociationReaction,
                type_1,
                type_2,
                delta_1,
                delta_2,
                min_state_1,
                max_state_1,
                min_state_2,
                max_state_2,
                cutoff,
                fpl,
                rate
            )
            self.cxxclass.set_reaction_cutoff(self, ReactionCutoffStaticLocal(cutoff))


    def add_postprocess(self, post_process, reactant_switch=0):
        """Add new post process to the reaction.
            Args:
                post_process: The post process object.
                reactant_switch: Which reactant (A, B or both) should be affected.
                    0 - both, 1 - reactant A, 2 - reactant B.
        """
        if pmi.workerIsActive():
            self.cxxclass.add_postprocess(self, post_process, reactant_switch)

if pmi.isController:
    class ChemicalReaction(Extension):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls='espressopp.integrator.ChemicalReactionLocal',
            pmiproperty=('interval','nearest_mode'),
            pmicall=(
                'add_reaction', 'get_timers'
                )
            )

    class PostProcessChangeProperty:
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls='espressopp.integrator.PostProcessChangePropertyLocal',
            pmicall=('add_change_property', 'remove_change_property')
        )

    class PostProcessChangeNeighboursProperty:
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls='espressopp.integrator.PostProcessChangeNeighboursPropertyLocal',
            pmicall=('add_change_property', )
        )

    class PostProcessChangePropertyOnState:
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls='espressopp.integrator.PostProcessChangePropertyOnStateLocal',
            pmicall=('add_change_property', )
        )

    class PostProcessRemoveBond:
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(cls='espressopp.integrator.PostProcessRemoveBondLocal')

    class PostProcessRemoveNeighbourBond:
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls='espressopp.integrator.PostProcessRemoveNeighbourBondLocal',
            pmicall=('add_bond_to_remove', )
        )

    class ReactionCutoffStatic:
        __metaclass__ = pmi.Proxy
        pmiproxydefs = {
            'cls': 'espressopp.integrator.ReactionCutoffStaticLocal',
            'pmiproperty': ('min_cutoff', 'max_cutoff')
        }

    class ReactionCutoffRandom:
        __metaclass__ = pmi.Proxy
        pmiproxydefs = {
            'cls': 'espressopp.integrator.ReactionCutoffRandomLocal',
            'pmiproperty': ('sigma', 'eq_distance')
        }

    class Reaction:
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls='espressopp.integrator.ReactionLocal',
            pmicall=(
                'add_postprocess',
                'set_reaction_cutoff',
                'get_reaction_cutoff'
            ),
            pmiproperty=(
                'type_1',
                'type_2',
                'delta_1',
                'delta_2',
                'min_state_1',
                'max_state_1',
                'min_state_2',
                'max_state_2',
                'intraresidual',
                'active',
                'cutoff',
                'rate',
                'is_virtual'
                )
            )

    class RestrictReaction:
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls='espressopp.integrator.RestrictReactionLocal',
            pmicall=(
                'add_postprocess',
                'set_reaction_cutoff',
                'get_reaction_cutoff',
                'define_connection'
            ),
            pmiproperty=(
                'type_1',
                'type_2',
                'delta_1',
                'delta_2',
                'min_state_1',
                'max_state_1',
                'min_state_2',
                'max_state_2',
                'intraresidual',
                'active',
                'cutoff',
                'rate',
                'revert',
            )
        )

    class DissociationReaction:
            __metaclass__ = pmi.Proxy
            pmiproxydefs = dict(
                cls='espressopp.integrator.DissociationReactionLocal',
                pmicall=(
                    'add_postprocess',
                ),
                pmiproperty=(
                    'type_1',
                    'type_2',
                    'delta_1',
                    'delta_2',
                    'min_state_1',
                    'max_state_1',
                    'min_state_2',
                    'max_state_2',
                    'cutoff',
                    'active',
                    'rate',
                    'diss_rate'
                    )
                )
