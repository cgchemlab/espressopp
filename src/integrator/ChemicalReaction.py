"""
 Copyright (C) 2014
   Pierre de Buyl, Jakub Krajniak
 Copyright (C) 2012,2013
   Max Planck Institute for Polymer Research
 Copyright (C) 2008,2009,2010,2011
   Max-Planck-Institute for Polymer Research & Fraunhofer SCAI

 This file is part of ESPResSo++.

 ESPResSo++ is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.

 ESPResSo++ is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""


from espresso.esutil import cxxinit
from espresso import pmi

from espresso.integrator.Extension import *  # NOQA
from _espresso import integrator_ChemicalReaction
from _espresso import integrator_Reaction
from _espresso import integrator_SynthesisReaction


class ChemicalReactionLocal(ExtensionLocal, integrator_ChemicalReaction):
    """Chemical Reaction scheme."""
    def __init__(self, system, vl, fpl, domdec):
        if (not (pmi._PMIComm and pmi._PMIComm.isActive()) or
                pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup()):
            cxxinit(self, integrator_ChemicalReaction, system, vl, fpl, domdec)

    def addReaction(self, reaction):
        if (not (pmi._PMIComm and pmi._PMIComm.isActive()) or
                pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup()):
            self.cxxclass.addReaction(self, reaction)

    def removeReaction(self, reaction_id):
        if (not (pmi._PMIComm and pmi._PMIComm.isActive()) or
                pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup()):
            self.cxxclass.removeReaction(self, reaction_id)


class SynthesisReactionLocal(integrator_SynthesisReaction, integrator_Reaction):
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
                'addReaction',
                'removeReaction'
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
