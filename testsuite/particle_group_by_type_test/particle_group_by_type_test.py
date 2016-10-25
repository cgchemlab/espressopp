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


class TestParticleGroupByTypeTest(ut.TestCase):
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
        self.particle_group = espressopp.ParticleGroupByType(self.system.storage, self.integrator)
        self.system.storage.decompose()
        self.integrator.run(1)

    def test_empty_group(self):
        self.assertEqual(self.particle_group.size(), 0)
        self.assertEqual(self.particle_group.get_particle_ids(), [[]])

    def test_single_type(self):
        self.particle_group.add_type_id(1)
        self.assertEqual(self.particle_group.size(), 0)
        self.integrator.run(1)
        self.assertEqual(self.particle_group.size(), 2)

    def test_multiple_type(self):
        self.particle_group.add_type_id(0)
        self.particle_group.add_type_id(1)
        self.assertEqual(self.particle_group.size(), 0)
        self.integrator.run(1)
        self.assertEqual(self.particle_group.size(), 3)

    def test_change_type(self):
        self.particle_group.add_type_id(0)
        self.particle_group.add_type_id(1)
        self.assertEqual(self.particle_group.size(), 0)
        self.integrator.run(1)
        self.assertEqual(self.particle_group.size(), 3)
        self.system.storage.modifyParticle(1, 'type', 4)
        self.system.storage.decompose()
        self.integrator.run(1)
        self.assertEqual(self.particle_group.size(), 2)

if __name__ == '__main__':
    ut.main()