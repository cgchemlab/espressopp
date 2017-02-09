#  Copyright (C) 2016,
#      Jakub Krajniak (jkrajniak at gmail.com)
#  Copyright (C) 2014
#      Pierre de Buyl
#  Copyright (C) 2012,2013
#      Max Planck Institute for Polymer Research
#  Copyright (C) 2008,2009,2010,2011
#      Max-Planck-Institute for Polymer Research & Fraunhofer SCAI
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
*****************************************
espressopp.interaction.LennardJones93Wall
*****************************************

This class defines a Lennard-Jones 9-3 SingleParticlePotential in the direction x.

.. math:: V(r) = \epsilon \left( \left(\frac{\sigma}{r}\right)^9 - \left(\frac{\sigma}{r}\right)^3 \right)

where :math:`r` is the distance from the lower or upper wall in the x
direction. :math:`V(r)=0` after a distance `sigmaCutoff`.

The parameters have to be defined for every species present in the system with
`setParams` and can be retrieved with `getParams`.

The direction of the wall is defined in the constructor.

Example:

    >>> LJ93 = espressopp.interaction.LennardJones93Wall(direction=0)
    >>> LJ93.setParams(0, 6., 1., wall_cutoff)
    >>> SPLJ93 = espressopp.interaction.SingleParticleLennardJones93Wall(system, LJ93)
    >>> system.addInteraction(SPLJ93)


.. function:: espressopp.interaction.LennardJones93Wall(direction)

      :param direction: The direction of the wall.
      :type direction: int


.. function:: espressopp.interaction.LennardJones93Wall.getParams(type_var)

		:param type_var: 
		:type type_var: 
		:rtype: 

.. function:: espressopp.interaction.LennardJones93Wall.setParams(type_var, epsilon, sigma, sigmaCutoff, r0)

		:param type_var: 
		:param epsilon: 
		:param sigma: 
		:param sigmaCutoff: 
		:param r0: 
		:type type_var: 
		:type epsilon: 
		:type sigma: 
		:type sigmaCutoff: 
		:type r0: 

.. function:: espressopp.interaction.SingleParticleLennardJones93Wall(system, potential)

		:param system: 
		:param potential: 
		:type system: 
		:type potential: 

.. function:: espressopp.interaction.SingleParticleLennardJones93Wall.setPotential(potential)

		:param potential: 
		:type potential: 
"""
from espressopp import pmi
from espressopp.esutil import *

from espressopp.interaction.SingleParticlePotential import *
from espressopp.interaction.Interaction import *
from _espressopp import interaction_LennardJones93Wall, interaction_SingleParticleLennardJones93Wall


class LennardJones93WallLocal(SingleParticlePotentialLocal, interaction_LennardJones93Wall):

    def __init__(self, direction=0):
        if pmi.workerIsActive():
            cxxinit(self, interaction_LennardJones93Wall, direction)

    def getParams(self, type_var):
        if pmi.workerIsActive():
            return self.cxxclass.getParams(self, type_var)

    def setParams(self, type_var, epsilon, sigma, sigmaCutoff, r0):
        if pmi.workerIsActive():
            self.cxxclass.setParams(self, type_var, epsilon, sigma, sigmaCutoff, r0)

class SingleParticleLennardJones93WallLocal(InteractionLocal, interaction_SingleParticleLennardJones93Wall):

    def __init__(self, system, potential):
        if pmi.workerIsActive():
            cxxinit(self, interaction_SingleParticleLennardJones93Wall, system, potential)

    def setPotential(self, potential):
        if pmi.workerIsActive():
            self.cxxclass.setPotential(self, potential)

if pmi.isController:
    class LennardJones93Wall(SingleParticlePotential):
        'The LennardJones93Wall potential.'
        pmiproxydefs = dict(
            cls = 'espressopp.interaction.LennardJones93WallLocal',
            pmicall = ['setParams', 'getParams']
            )

    class SingleParticleLennardJones93Wall(Interaction):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls = 'espressopp.interaction.SingleParticleLennardJones93WallLocal',
            pmicall = ['setPotential']
            )
