#  Copyright (C) 2014-2015
#    Jakub Krajniak (jkrajniak@gmail.com)
#  Copyright (C) 2012,2013
#    Max Planck Institute for Polymer Research
#  Copyright (C) 2008,2009,2010,2011
#    Max-Planck-Institute for Polymer Research & Fraunhofer SCAI
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
**********************
**Chemical reactions**
**********************

This extension enables the rate-controlled stochastic curing of polymer
systems, either for chain growth of step growth, depending on the
parameters.

The variables typeA, typeB, minStateA, minStateB, maxStateA, maxStateB
control the particles that enter the curing reaction

.. math::
  A^a + B^b \\rightarrow A^{a+deltaA}-B^{b+deltaB}

where A and B may possess additional bonds not shown.

An extra bond is added between A and B whenever the state of A and B falls
into the defined range by variables min/max state.
The condition is as follow:

.. math::

   a >= minStateA \\land stateA < maxStateA

the same holds for the particle B. Both condition should match.
In addition if the intramolecular property is set to true (by default) then
the reaction only happen between heterogeneous molecules.

The reaction proceeds by testing for all possible (A,B) pairs and
selects them only at a given rate. It works in parallel, by gathering
first the successful pairs between neighboring CPUs and ensuring that
each particle enters only in one new bond per reaction step.

Example
#######

**Creating the integrator extension**

>>> ar = espressopp.integrator.ChemicalReaction(
>>>     system,
>>>     verletList,
>>>     fpl_a_b,
>>>     system.storage,
>>>     interval)

**Creates synthesis reaction**

>>> r_type_1 = espressopp.integrator.Reaction(
>>>     type_a=ar_type_M,
>>>     type_b=ar_type_B,
>>>     delta_a=1,
>>>     delta_b=1,
>>>     min_state_a=0,
>>>     min_state_b=0,
>>>     max_state_a=2,
>>>     max_state_b=1,
>>>     rate=1000,
>>>     cutoff=ar_cutoff
>>>     )

Add the reaction to the integrator extension

>>> ar.add_reaction(r_type_1)

Add the extension to the integrator

>>> integrator.addExtension(ar)

"""


from espressopp.esutil import cxxinit
from espressopp import pmi

from espressopp.integrator.Extension import *  # NOQA
from _espressopp import integrator_ChemicalReaction
from _espressopp import integrator_Reaction

from _espressopp import integrator_PostProcess
from _espressopp import integrator_PostProcessChangeProperty


class ChemicalReactionLocal(ExtensionLocal, integrator_ChemicalReaction):
    """Chemical Reaction integrator extension."""
    def __init__(self, system, vl, fpl, domdec, interval=None):
        """Chemical reaction extension.

        Args:
          system: The espressopp.System object.
          vl: The verlet list.
          fpl: The fixed pair list that will hold new bonds.
          domdec: The domain decomposition object.
          interval: The timestep where the extension will be run.
        """
        if (not (pmi._PMIComm and pmi._PMIComm.isActive()) or
                pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup()):
            cxxinit(self, integrator_ChemicalReaction, system, vl, fpl, domdec)
            if interval is not None:
                self.interval = interval

    def add_reaction(self, reaction):
        """Adds the reaction to the list."""
        if (not (pmi._PMIComm and pmi._PMIComm.isActive()) or
                pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup()):
            self.cxxclass.addReaction(self, reaction)

    def remove_reaction(self, reaction_id):
        """Removes the reaction from the list.

        Args:
          reaction_id: The id of the reaction to remove.
        """
        if (not (pmi._PMIComm and pmi._PMIComm.isActive()) or
                pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup()):
            self.cxxclass.removeReaction(self, reaction_id)


class PostProcessChangePropertyLocal(integrator_PostProcessChangeProperty,
                                     integrator_PostProcess):
    """Post process of reaction that changes particle property."""
    def __init__(self):
        if (not (pmi._PMIComm and pmi._PMIComm.isActive()) or
                pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup()):
            cxxinit(self, integrator_PostProcessChangeProperty)

    def add_change_property(self, type_id, prop):
        if (not (pmi._PMIComm and pmi._PMIComm.isActive()) or
                pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup()):
            self.cxxclass.add_change_property(self, type_id, prop)

    def remove_change_property(self, type_id):
        if (not (pmi._PMIComm and pmi._PMIComm.isActive()) or
                pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup()):
            self.cxxclass.remove_change_property(self, type_id)


class ReactionLocal(integrator_Reaction):
    """Synthesis reaction."""
    def __init__(self, type_1, type_2, delta_1, delta_2, min_state_1, max_state_1,
                 min_state_2, max_state_2, cutoff, rate, intramolecular=False):
        if (not (pmi._PMIComm and pmi._PMIComm.isActive()) or
                pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup()):
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
                intramolecular
            )

    def add_postprocess(self, post_process):
        if (not (pmi._PMIComm and pmi._PMIComm.isActive()) or
                pmi._MPIcomm.rank in pmi._PMIComm.GetMPIcpugroup()):
            self.cxxclass.add_postprocess(self, post_process)


if pmi.isController:
    class ChemicalReaction(Extension):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls='espressopp.integrator.ChemicalReactionLocal',
            pmiproperty=('interval',),
            pmicall=(
                'add_reaction',
                'remove_reaction'
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

    class PostProcessUpdateExcludeList:
        __metaclass__ = pmi.Proxy
        pmiproxydefs = {
            'cls': 'espressopp.integrator.PostProcessUpdateExcludeListLocal'
        }

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
                'intramolecular'
                )
            )
