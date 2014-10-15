#  Copyright (C) 2014
#    Pierre de Buyl
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

>>> ar = espresso.integrator.ChemicalReaction(
>>>     system,
>>>     verletList,
>>>     fpl_a_b,
>>>     system.storage,
>>>     interval)

**Creates synthesis reaction**

>>> r_type_1 = espresso.integrator.SynthesisReaction(
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


from espresso.esutil import cxxinit
from espresso import pmi

from espresso.integrator.Extension import *  # NOQA
from _espresso import integrator_ChemicalReaction
from _espresso import integrator_Reaction
from _espresso import integrator_SynthesisReaction


class ChemicalReactionLocal(ExtensionLocal, integrator_ChemicalReaction):
    """Chemical Reaction integrator extension."""
    def __init__(self, system, vl, fpl, domdec, interval=None):
        """Chemical reaction extension.
        
        Args:
          system: The espresso.System object.
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


class SynthesisReactionLocal(integrator_SynthesisReaction, integrator_Reaction):
    """Synthesis reaction."""
    def __init__(self, type_a, type_b, delta_a, delta_b, min_state_a, min_state_b,
                 max_state_a, max_state_b, cutoff, rate, intramolecular=False):
        if (not (pmi._PMIComm and pmi._PMIComm.isActive()) or
                pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup()):
            cxxinit(
                self,
                integrator_SynthesisReaction,
                type_a,
                type_b,
                delta_a,
                delta_b,
                min_state_a,
                min_state_b,
                max_state_a,
                max_state_b,
                cutoff,
                rate,
                intramolecular
            )


if pmi.isController:
    class ChemicalReaction(Extension):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls='espresso.integrator.ChemicalReactionLocal',
            pmiproperty=('interval',),
            pmicall=(
                'add_reaction',
                'remove_reaction'
                )
            )

    class SynthesisReaction:
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls='espresso.integrator.SynthesisReactionLocal',
            pmiproperty=(
                'type_a',
                'type_b',
                'delta_a',
                'delta_b',
                'min_state_a',
                'max_state_a',
                'min_state_b',
                'max_state_b',
                'rate',
                'cutoff',
                'intramolecular'
                )
            )
