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
**espressopp.integrator.ChemicalReaction**
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

"""


from espressopp.esutil import cxxinit
from espressopp import pmi

from espressopp.integrator.Extension import *  # NOQA
from _espressopp import integrator_ChemicalReaction
from _espressopp import integrator_Reaction
from _espressopp import integrator_DissociationReaction

from _espressopp import integrator_ChemicalReactionPostProcess
from _espressopp import integrator_PostProcessChangeProperty
from _espressopp import integrator_PostProcessRemoveBond


class ChemicalReactionLocal(ExtensionLocal, integrator_ChemicalReaction):
    """Chemical Reaction integrator extension."""
    def __init__(self, system, vl, domdec, interval=None):
        """Chemical reaction extension.

        Args:
          system: The espressopp.System object.
          vl: The verlet list.
          domdec: The domain decomposition object.
          interval: The timestep where the extension will be run.
        """
        if pmi.workerIsActive():
            cxxinit(self, integrator_ChemicalReaction, system, vl, domdec)
            if interval is not None:
                self.interval = interval

    def add_reaction(self, reaction):
        """Adds the reaction to the list."""
        if pmi.workerIsActive():
            self.cxxclass.add_reaction(self, reaction)


class PostProcessChangePropertyLocal(integrator_PostProcessChangeProperty,
                                     integrator_ChemicalReactionPostProcess):
    """Post process of reaction that changes particle property."""
    def __init__(self):
        if pmi.workerIsActive():
            cxxinit(self, integrator_PostProcessChangeProperty)

    def add_change_property(self, type_id, prop):
        if pmi.workerIsActive():
            self.cxxclass.add_change_property(self, type_id, prop)

    def remove_change_property(self, type_id):
        if pmi.workerIsActive():
            self.cxxclass.remove_change_property(self, type_id)

class PostProcessRemoveBondLocal(integrator_PostProcessRemoveBond, integrator_ChemicalReactionPostProcess):
    """Post process of reaction, remove some bonds."""
    def __init__(self, fpl, number_of_bonds):
        if pmi.workerIsActive():
            cxxinit(self, integrator_PostProcessRemoveBond, fpl, number_of_bonds)


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
                cutoff,
                rate,
                fpl,
                False
            )

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
                rate,
                fpl
            )

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
            pmiproperty=('interval',),
            pmicall=(
                'add_reaction',
                )
            )

    class PostProcessChangeProperty:
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls='espressopp.integrator.PostProcessChangePropertyLocal',
            pmicall=('add_change_property', 'remove_change_property')
        )

    class PostProcessUpdateResId:
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls='espressopp.integrator.PostProcessUpdateResIdLocal',
            pmicall=('add_molecule_size',)
        )

    class PostProcessRemoveBond:
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(cls='espressopp.integrator.PostProcessRemoveBondLocal')

    class Reaction:
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls='espressopp.integrator.ReactionLocal',
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
                'rate',
                'cutoff',
                'min_cutoff',
                'intramolecular',
                'active'
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
                    'rate',
                    'cutoff',
                    'min_cutoff',
                    'diss_rate',
                    'active'
                    )
                )