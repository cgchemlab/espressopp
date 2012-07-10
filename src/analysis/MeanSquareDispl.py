from espresso.esutil import cxxinit
from espresso import pmi

from espresso.analysis.ConfigsParticleDecomp import *
from _espresso import analysis_MeanSquareDispl

class MeanSquareDisplLocal(ConfigsParticleDecompLocal, analysis_MeanSquareDispl):
    'The (local) compute autocorrelation f.'
    def __init__(self, system):
      cxxinit(self, analysis_MeanSquareDispl, system)
      
if pmi.isController:
  class MeanSquareDispl(ConfigsParticleDecomp):
    __metaclass__ = pmi.Proxy
    pmiproxydefs = dict(
      cls =  'espresso.analysis.MeanSquareDisplLocal'
    )
