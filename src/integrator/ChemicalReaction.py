from espresso.esutil import cxxinit
from espresso import pmi

from espresso.integrator.Extension import *  #NOQA
from _espresso import integrator_ChemicalReaction


class ChemicalReactionLocal(ExtensionLocal, integrator_ChemicalReaction):
    """Chemical Reaction scheme."""
    def __init__(self, system, vl, fpl, domdec):
        if (not (pmi._PMIComm and pmi._PMIComm.isActive()) or
                pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup()):
            cxxinit(self, integrator_ChemicalReaction, system, vl, fpl, domdec)


class ReactionLocal():
    """Base class for Reaction scheme."""
    def __init__(self):
        if (not (pmi._PMIComm and pmi._PMIComm.isActive()) or
                pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup()):
            cxxinit(self)


class SynthesisReactionLocal(ReactionLocal):
    """Synthesis reaction scheme."""
    def __init__(self):
        if (not (pmi._PMIComm and pmi._PMIComm.isActive()) or
                pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup()):
            cxxinit(self)


if pmi.isController:
    class ChemicalReaction(Extension):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls='espresso.integrator.ChemicalReactionLocal',
            pmiproperty=['interval']
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
                'cutoff'
                ]
            )

    class SynthesisReaction(Reaction):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(cls='espresso.integrator.SynthesisReactionLocal')
