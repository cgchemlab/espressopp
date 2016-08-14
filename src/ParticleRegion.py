#  Copyright (C) 2016
#      Jakub Krajniak (c) (jkrajniak at gmail.com)
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
*****************************
**espressopp.ParticleRegion**
*****************************

Defines a geometrical block and gets the list of particles in that region, updated
dynamicly during the simulation. There is no limit on number of such regions in the
simulation box.
The **ParticleRegion** can store particles of any type or only of the type
defined byt the **add_type_id**.

The particle is in the region whenever :math:`p > l ^ p < r` where
:math:`p` is a position of particle, :math:`l`: is the position of left-bottom-front vertex
of the block and :math:`r` is the position of right-top-rear vertex.

**ParticleRegion** inherits interface from **ParticleGroup**, therefore whenever
**ParticleGroup** is accepted then **ParticleRegion** can also be passed.

The region can be extended or shrinked during the simulation by defining the speed of it. The update
is always made on the right top bottom coordinate of the region and the velocity can be positive or negative. By
default it is 0.0 in all directions.

Example of usage
==================

>>> pr = espressopp.ParticleRegion(storage, espressopp.Real3D(-1, -1, 5), espressopp.Real3D(-1, -1, 11))


.. class:: espressopp.ParticleRegion

      The particle region class.

      .. method:: espressopp.ParticleRegion(storage, left, right)

            Define the particle region.

            :param storage: The storage object.
		    :type storage: espressopp.storage.Storage
		    :param left: The left-bottom-front vertex.
		    :type left: espressopp.Real3D
		    :param right: The right-top-rear vertex.
		    :type right: espressopp.Real3D

	  .. method:: espressopp.ParticleRegion.set_v(vx, vy, vz)

            Set the velocity of which the region will change the dimension. The velocity
            can be positive or negative.

            :param vx: The velocity in x-direction.
            :type vx: float
            :param vy: The velocity in y-direction.
            :type vy: float
            :param vz: The velocity in z-direction.
            :type vz: float

	  .. method:: espressopp.ParticleRegion.get_v()

            Returns the velocity of the changes in region.

            :rtype: tuple

      .. method:: espressopp.ParticleRegion.has(pid)

            Check if given particle of *pid* is in the region.

            :param pid:
            :type pid:
            :rtype: bool

      .. method:: espressopp.ParticleRegion.show()

            Debug purpose, print all particles in the region.


      .. method:: espressopp.ParticleRegion.size()

            Gets the number of particles in the region.

            :rtype: int


      .. method:: espressopp.ParticleRegion.add_type_id(type_id)

            Stores only particles of given type.

            :param type_id: Type id to add
            :type type_id: int


      .. method:: espressopp.ParticleRegion.remove_type_id(type_id)

            Removes type from the list of observed types.

            :param type_id: Type id to remove.
            :type type_id: int


.. function:: espressopp.ParticleRegion(storage, left, right)

		:param storage: The storage object.
		:type storage: espressopp.storage.Storage
		:param left: The left-bottom-front vertex.
		:type left: espressopp.Real3D
		:param right: The right-top-rear vertex.
		:type right: espressopp.Real3D

.. function:: espressopp.ParticleRegion.has(pid)
   Check if given particle of *pid* is in the region.

		:param pid: 
		:type pid: 
		:rtype: 

.. function:: espressopp.ParticleRegion.show()
   Debug purpose, print all particles in the region.
		:rtype: 

.. function:: espressopp.ParticleRegion.size()
   Gets the number of particles in the region.
		:rtype:


.. function:: espressopp.ParticleRegion.add_type_id(type_id)
   Stores only particles of given type.

        :param type_id: Type id to add
        :type type_id: int


.. function:: espressopp.ParticleRegion.remove_type_id(type_id)
   Removes type from the list of observed types.

        :param type_id: Type id to remove.
        :type type_id: int

"""
import _espressopp
from espressopp.esutil import cxxinit

import pmi


class ParticleRegionLocal(_espressopp.ParticleRegion):
    def __init__(self, storage, integrator, left_bottom, right_top):
        if pmi.workerIsActive():
            cxxinit(self, _espressopp.ParticleRegion, storage, integrator)
            self.cxxclass.define_region(self, left_bottom, right_top)

    def add_type_id(self, type_id):
        if pmi.workerIsActive():
            self.cxxclass.add_type_id(self, type_id)

    def remove_type_id(self, type_id):
        if pmi.workerIsActive():
            self.cxxclass.remove_type_id(self, type_id)

    def show(self):
        if pmi.workerIsActive():
            self.cxxclass.show(self)

    def has(self, pid):
        if pmi.workerIsActive():
            return self.cxxclass.has(self, pid)

    def size(self):
        if pmi.workerIsActive():
            return self.cxxclass.size(self)

    def set_v(self, vx, vy, vz, left_right='left'):
        if pmi.workerIsActive():
            self.cxxclass.set_v(self, vx, vy, vz, left_right=='left')

    def get_v(self, left_right='left'):
        if pmi.workerIsActive():
            return self.cxxclass.get_v(self, left_right=='left')

    def get_region(self):
        if pmi.workerIsActive():
            return self.cxxclass.get_region(self)


if pmi.isController:
    class ParticleRegion(object):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls='espressopp.ParticleRegionLocal',
            pmiinvoke=['get_particle_ids'],
            pmicall=['show', 'has', 'size', 'define_region', 'add_type_id', 'remove_type_id', 'set_v', 'get_v',
                     'get_region']
        )
