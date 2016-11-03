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


r"""
*****************************************************
**espressopp.interaction.LennardJonesSoftCoreLambda**
*****************************************************
.. math::
	V(r) = 4 \varepsilon \left[ \left( \frac{\sigma}{r_{eff}} \right)^{12} -
	\left( \frac{\sigma}{r_{eff}} \right)^{6} \right]

.. math::
    r_{eff} = \left[ \alpha \sigma^6 (1-w(r))^P + r^6 \right]^{1/6}

Reference:
Peters, J. H. et al. Phys. Rev. E - Stat. Nonlinear, Soft Matter Phys. 94, 1-19 (2016).

"""
from espressopp import pmi, infinity
from espressopp.esutil import *

from espressopp.interaction.Potential import *
from espressopp.interaction.Interaction import *
from _espressopp import interaction_LennardJonesSoftCoreLambda, \
                        interaction_VerletListLennardJonesSoftCoreLambda

class LennardJonesSoftCoreLambdaLocal(PotentialLocal, interaction_LennardJonesSoftCoreLambda):
    def __init__(self, epsilon=1.0, sigma=1.0, cutoff=infinity, shift="auto", alpha=1.0, weight_power=1.0):
        """Initialize the local Lennard Jones object."""
        if pmi.workerIsActive():
            if shift =="auto":
                cxxinit(self, interaction_LennardJonesSoftCoreLambda,
                        epsilon, sigma, cutoff, alpha, weight_power)
            else:
                cxxinit(self, interaction_LennardJonesSoftCoreLambda,
                        epsilon, sigma, cutoff, shift, alpha, weight_power)


class VerletListLennardJonesSoftCoreLambdaLocal(InteractionLocal, interaction_VerletListLennardJonesSoftCoreLambda):
    def __init__(self, vl):
        if pmi.workerIsActive():
            cxxinit(self, interaction_VerletListLennardJonesSoftCoreLambda, vl)

    def setPotential(self, type1, type2, potential):
        if pmi.workerIsActive():
            self.cxxclass.setPotential(self, type1, type2, potential)

    def getPotential(self, type1, type2):
        if pmi.workerIsActive():
            return self.cxxclass.getPotential(self, type1, type2)

    def getVerletListLocal(self):
        if pmi.workerIsActive():
            return self.cxxclass.getVerletList(self)


if pmi.isController:
    class LennardJonesSoftCoreLambda(Potential):
        'The Lennard-Jones potential.'
        pmiproxydefs = dict(
            cls = 'espressopp.interaction.LennardJonesSoftCoreLambdaLocal',
            pmiproperty = ['epsilon', 'sigma', 'alpha', 'weight_power']
            )

    class VerletListLennardJonesSoftCoreLambda(Interaction):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls =  'espressopp.interaction.VerletListLennardJonesSoftCoreLambdaLocal',
            pmicall = ['setPotential', 'getPotential', 'getVerletList']
            )
