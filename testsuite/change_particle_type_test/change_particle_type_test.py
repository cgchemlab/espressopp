#!/usr/bin/env python
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

import espressopp
import unittest as ut


class TestChangeParticleTypeTest(ut.TestCase):
    def setUp(self):
        # super(TestParticleGroupByTypeTest, self).setUp()
        self.system, self.integrator = espressopp.standard_system.Minimal(
            0, (10., 10., 10.), dt=0.01)
        self.particle_list = [
            (1, 0, espressopp.Real3D(1, 1, 1)),
            (2, 1, espressopp.Real3D(1, 1, 6)),
            (3, 1, espressopp.Real3D(1, 1, 4.7)),
            (4, 2, espressopp.Real3D(1, 1, 10))]
        self.part_prop = ['id', 'type', 'pos']
        self.system.storage.addParticles(self.particle_list, *self.part_prop)
        self.system.storage.decompose()

    def test_interval(self):
        change_type = espressopp.integrator.ChangeParticleType(
            self.system, 10, 2, 1, 3)
        self.integrator.addExtension(change_type)
        self.integrator.run(5)
        self.assertEqual(self.system.storage.getParticle(2).type, 1)

    def test_change_type(self):
        change_type = espressopp.integrator.ChangeParticleType(
            self.system, 10, 2, 1, 3)
        self.integrator.addExtension(change_type)
        self.integrator.run(15)
        self.assertEqual(self.system.storage.getParticle(1).type, 0)
        self.assertEqual(self.system.storage.getParticle(2).type, 3)  # change type: 1->3
        self.assertEqual(self.system.storage.getParticle(3).type, 3)  # change type: 1->3
        self.assertEqual(self.system.storage.getParticle(4).type, 2)

if __name__ == '__main__':
    ut.main()