from espresso.esutil import cxxinit
from espresso import pmi

from espresso.integrator.Extension import *  #NOQA
from _espresso import integrator_ChemicalReaction
from _espresso import integrator_Reaction


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


class ReactionLocal(integrator_Reaction):
    """Base class for Reaction scheme."""
    def __init__(self, type_a, type_b, delta_a, delta_b, min_state_a, min_state_b,
                 max_state_a, max_state_b, cutoff, rate, intramolecular=False):
        if (not (pmi._PMIComm and pmi._PMIComm.isActive()) or
                pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup()):
            cxxinit(
                self,
                integrator_Reaction,
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

"""
class SynthesisReactionLocal(integrator_SynthesisReaction, ReactionLocal):
    def __init__(self):
        if (not (pmi._PMIComm and pmi._PMIComm.isActive()) or
                 pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup()):
            cxxinit(self, integrator_SynthesisReaction)
"""

if pmi.isController:
    class ChemicalReaction(Extension):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls='espresso.integrator.ChemicalReactionLocal',
            pmiproperty=['interval'],
            pmicall=['addReaction', 'removeReaction']
            )

    class Reaction():
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls='espresso.integrator.ReactionLocal',
            pmiproperty=[
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
                ]
            )
        """
        class SynthesisReaction():
            __metaclass__ = pmi.Proxy
            pmiproxydefs = dict(
                cls='espresso.integrator.SynthesisReactionLocal',
                )
        """
