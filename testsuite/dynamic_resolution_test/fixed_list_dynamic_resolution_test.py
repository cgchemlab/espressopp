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


import unittest as ut
import espressopp


class ESPPTestCase(ut.TestCase):
    def setUp(self):
        self.system, self.integrator = espressopp.standard_system.Minimal(
            0, (10., 10., 10.))
        particles = [
            (1, 0, espressopp.Real3D(5, 5, 5)),
            (2, 0, espressopp.Real3D(5, 6.97, 5)),
            (3, 1, espressopp.Real3D(5, 8.94, 5)),
            (4, 1, espressopp.Real3D(5, 10.94, 5))
        ]
        self.part_prop = ('id', 'type', 'pos')
        self.system.storage.addParticles(particles, *self.part_prop)
        self.fixed_dynamic_resolution = espressopp.integrator.FixedListDynamicResolution(self.system)
        self.integrator.addExtension(self.fixed_dynamic_resolution)


class TestFixedPairListLambda(ESPPTestCase):
    def test_harmonic_energy(self):
        fpl = espressopp.FixedPairListLambda(self.system.storage)
        fpl.addBonds([(1, 2)])
        harmonic = espressopp.interaction.Harmonic(K=30, r0=0.97)
        interaction = espressopp.interaction.FixedPairListLambdaHarmonic(self.system, fpl, harmonic)
        self.assertAlmostEqual(interaction.computeEnergy(), 30.0)
        self.fixed_dynamic_resolution.register_pair_list(fpl, -0.5)
        self.integrator.run(1)
        self.assertAlmostEqual(interaction.computeEnergy(), 15.0)
        self.fixed_dynamic_resolution.update_pair_list(0, 0.5)  # change rate of fpl, index 0
        self.integrator.run(1)
        self.assertAlmostEqual(interaction.computeEnergy(), 30.0)


class TestFixedTripleListLambda(ESPPTestCase):
    def test_angular_harmonic(self):
        ftl_lambda = espressopp.FixedTripleListLambda(self.system.storage)
        ftl_lambda.addTriples([(1, 2, 3)])
        ftl = espressopp.FixedTripleList(self.system.storage)
        ftl.addTriples([(1, 2, 3)])
        harmonic = espressopp.interaction.AngularHarmonic()
        interactionLambda = espressopp.interaction.FixedTripleListLambdaAngularHarmonic(self.system, ftl_lambda, harmonic)
        interaction = espressopp.interaction.FixedTripleListAngularHarmonic(self.system, ftl, harmonic)
        self.assertAlmostEqual(interaction.computeEnergy(), interactionLambda.computeEnergy())

    def test_angular_harmonic_lambda(self):
        ftl_lambda = espressopp.FixedTripleListLambda(self.system.storage)
        ftl_lambda.addTriples([(1, 2, 3)])
        ftl = espressopp.FixedTripleList(self.system.storage)
        ftl.addTriples([(1, 2, 3)])
        harmonic = espressopp.interaction.AngularHarmonic()
        interactionLambda = espressopp.interaction.FixedTripleListLambdaAngularHarmonic(self.system, ftl_lambda, harmonic)
        interaction = espressopp.interaction.FixedTripleListAngularHarmonic(self.system, ftl, harmonic)

        self.fixed_dynamic_resolution.register_triple_list(ftl_lambda, -0.5)
        self.integrator.run(1)
        self.assertAlmostEqual(0.5*interaction.computeEnergy(), interactionLambda.computeEnergy())


class TestFixedQuadrupleListLambda(ESPPTestCase):
    def setUp(self):
        super(TestFixedQuadrupleListLambda, self).setUp()
        for i in range(1, 5):
            self.system.storage.modifyParticle(i, 'pos', self.system.bc.getRandomPos())
        self.potential = espressopp.interaction.DihedralRB(1.0, 0.5, 0.25, 0.125)

    def test_dihedral_harmonic(self):
        ftl_lambda = espressopp.FixedQuadrupleListLambda(self.system.storage)
        ftl_lambda.addQuadruples([(1, 2, 3, 4)])
        ftl = espressopp.FixedQuadrupleList(self.system.storage)
        ftl.addQuadruples([(1, 2, 3, 4)])

        interactionLambda = espressopp.interaction.FixedQuadrupleListLambdaDihedralRB(self.system, ftl_lambda, self.potential)
        interaction = espressopp.interaction.FixedQuadrupleListDihedralRB(self.system, ftl, self.potential)
        self.assertAlmostEqual(interaction.computeEnergy(), interactionLambda.computeEnergy())

    def test_dihedral_harmonic_lambda(self):
        ftl_lambda = espressopp.FixedQuadrupleListLambda(self.system.storage)
        ftl_lambda.addQuadruples([(1, 2, 3, 4), (2, 1, 3, 4)])
        ftl = espressopp.FixedQuadrupleList(self.system.storage)
        ftl.addQuadruples([(1, 2, 3, 4), (2, 1, 3, 4)])
        interactionLambda = espressopp.interaction.FixedQuadrupleListLambdaDihedralRB(
            self.system, ftl_lambda, self.potential)
        interaction = espressopp.interaction.FixedQuadrupleListDihedralRB(self.system, ftl, self.potential)

        self.fixed_dynamic_resolution.register_quadruple_list(ftl_lambda, -0.5)
        self.integrator.run(1)
        self.assertAlmostEqual(0.5*interaction.computeEnergy(), interactionLambda.computeEnergy())


if __name__ == '__main__':
    ut.main()