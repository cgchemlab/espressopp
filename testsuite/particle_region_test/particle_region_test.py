#!/usr/bin/env python
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

import espressopp
import unittest as ut


class TestFixedPairListTypesTabulated(ut.TestCase):
    def setUp(self):
        self.system, self.integrator = espressopp.standard_system.Minimal(
            0, (10., 10., 10.))
        self.system.storage.addParticle(1, espressopp.Real3D(1, 1, 1))
        self.system.storage.addParticle(2, espressopp.Real3D(1, 1, 6))
        self.system.storage.addParticle(3, espressopp.Real3D(1, 1, 4.7))
        self.system.storage.modifyParticle(3, 'v', espressopp.Real3D(0, 0, 0.5))
        self.system.storage.addParticle(4, espressopp.Real3D(1, 1, 10))
        self.particle_region = espressopp.ParticleRegion(
            self.system.storage,
            espressopp.Real3D(-1, -1, 5),
            espressopp.Real3D(11, 11, 11))
        self.system.storage.decompose()

    def test_static_region(self):
        self.assertEqual(self.particle_region.size(), 2)
        self.assertEqual(self.particle_region.get_particle_ids(), [[2, 4]])

    def test_particle_gets_into_region(self):
        # Now run integrator, particle 3 will move into region after while
        for i in range(10):
            self.integrator.run(100)

        self.assertEqual(self.particle_region.size(), 3)  # Particle 3 in the region
        self.assertEqual(self.particle_region.get_particle_ids(), [[2, 3, 4]])


if __name__ == '__main__':
    ut.main()