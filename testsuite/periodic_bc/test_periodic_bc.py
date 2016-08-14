#!/usr/bin/env python

import espressopp
import mpi4py.MPI as MPI

import unittest


class TestPeriodicBC(unittest.TestCase):
    def setUp(self):
        system = espressopp.System()
        system.skin = 0.3
        system.comm = MPI.COMM_WORLD
        self.system = system
        self.box = (10.0, 10.0, 10.0)
        NCPUs = espressopp.MPI.COMM_WORLD.size
        # calculate a regular 3D grid according to the number of CPUs available
        self.nodeGrid = espressopp.tools.decomp.nodeGrid(NCPUs)
        # calculate a 3D subgrid to speed up verlet list builds and communication
        self.cellGrid = espressopp.tools.decomp.cellGrid(self.box, self.nodeGrid, 3.0, system.skin)

    def test_periodic_system(self):
        """Full periodic system, test the distance method."""
        self.system.bc = espressopp.bc.FreeOrthorhombicBC(self.system.rng, (10, 10, 10))
        self.system.storage = espressopp.storage.DomainDecompositionFree(self.system, self.nodeGrid, self.cellGrid)
        pos_diff = self.system.bc.getMinimumImageVector(
            espressopp.Real3D(9.5, 9.0, 2.0),
            espressopp.Real3D(0.5, 2.0, 2.0))
        self.assertEqual(pos_diff, (-1.0, -3.0, 0.0))

    def test_non_periodic_system(self):
        """Switch off periodicity in x-direction"""
        self.system.bc = espressopp.bc.FreeOrthorhombicBC(
            self.system.rng, (10, 10, 10),
            [False, True, True])
        self.system.storage = espressopp.storage.DomainDecompositionFree(self.system, self.nodeGrid, self.cellGrid)
        pos_diff = self.system.bc.getMinimumImageVector(
            espressopp.Real3D(9.5, 9.0, 2.0),
            espressopp.Real3D(0.5, 2.0, 2.0))
        self.assertEqual(pos_diff, (9.0, -3.0, 0.0))

    def test_potential_periodic_bounduary(self):
        self.system.bc = espressopp.bc.FreeOrthorhombicBC(
            self.system.rng, (10, 10, 10), [True, True, True])
        self.system.storage = espressopp.storage.DomainDecompositionFree(self.system, self.nodeGrid, self.cellGrid)
        particle_list = [
            (1, espressopp.Real3D(9.5, 2.0, 2.0), 1.0),
            (2, espressopp.Real3D(0.5, 2.0, 2.0), 1.0),
        ]
        self.system.storage.addParticles(particle_list, 'id', 'pos', 'mass')
        self.system.storage.decompose()

        vl = espressopp.VerletList(self.system, cutoff=2.5)
        lj = espressopp.interaction.LennardJones(sigma=1.0, epsilon=1.0, cutoff=2.5)
        interaction = espressopp.interaction.VerletListLennardJones(vl)
        interaction.setPotential(type1=0, type2=0, potential=lj)

        self.system.addInteraction(interaction)
        interaction_energy = interaction.computeEnergy()
        self.assertEqual(interaction_energy, 0.016316891136000003)

    def test_potential_non_periodic_bounduary(self):
        self.system.bc = espressopp.bc.FreeOrthorhombicBC(
            self.system.rng, (10, 10, 10), [False, True, True])
        self.system.storage = espressopp.storage.DomainDecompositionFree(self.system, self.nodeGrid, self.cellGrid)
        particle_list = [
            (1, espressopp.Real3D(9.5, 2.0, 2.0), 1.0),
            (2, espressopp.Real3D(0.5, 2.0, 2.0), 1.0),
        ]
        self.system.storage.addParticles(particle_list, 'id', 'pos', 'mass')
        self.system.storage.decompose()

        vl = espressopp.VerletList(self.system, cutoff=2.5)
        lj = espressopp.interaction.LennardJones(sigma=1.0, epsilon=1.0, cutoff=2.5)
        interaction = espressopp.interaction.VerletListLennardJones(vl)
        interaction.setPotential(type1=0, type2=0, potential=lj)

        self.system.addInteraction(interaction)
        interaction_energy = interaction.computeEnergy()
        # Particles does not interact across the box in x-direction
        self.assertEqual(interaction_energy, 0.0)

    def test_potential_non_periodic_inner(self):
        #espressopp.PLogger.set('DomainDecompositionFree', 'TRACE')
        self.system.bc = espressopp.bc.FreeOrthorhombicBC(
            self.system.rng, (10, 10, 10), [False, False, False])
        self.system.storage = espressopp.storage.DomainDecompositionFree(self.system, self.nodeGrid, self.cellGrid)
        particle_list = [
            (1, espressopp.Real3D(4.0, 2.0, 2.0), 1.0),
            (2, espressopp.Real3D(5.5, 2.0, 2.0), 1.0),
        ]
        self.system.storage.addParticles(particle_list, 'id', 'pos', 'mass')
        self.system.storage.decompose()

        vl = espressopp.VerletList(self.system, cutoff=2.5)
        lj = espressopp.interaction.LennardJones(sigma=1.0, epsilon=1.0, cutoff=2.5)
        interaction = espressopp.interaction.VerletListLennardJones(vl)
        interaction.setPotential(type1=0, type2=0, potential=lj)

        self.system.addInteraction(interaction)
        interaction_energy = interaction.computeEnergy()

        # Particles does not interact across the box in x-direction
        self.assertEqual(interaction_energy, -0.30401970314257465)

    def test_potential_periodic_inner(self):
        #espressopp.PLogger.set('DomainDecompositionFree', 'TRACE')
        self.system.bc = espressopp.bc.FreeOrthorhombicBC(
            self.system.rng, (10, 10, 10), [True, True, True])
        self.system.storage = espressopp.storage.DomainDecompositionFree(self.system, self.nodeGrid, self.cellGrid)
        particle_list = [
            (1, espressopp.Real3D(4.0, 2.0, 2.0), 1.0),
            (2, espressopp.Real3D(5.5, 2.0, 2.0), 1.0),
        ]
        self.system.storage.addParticles(particle_list, 'id', 'pos', 'mass')
        self.system.storage.decompose()

        vl = espressopp.VerletList(self.system, cutoff=2.5)
        lj = espressopp.interaction.LennardJones(sigma=1.0, epsilon=1.0, cutoff=2.5)
        interaction = espressopp.interaction.VerletListLennardJones(vl)
        interaction.setPotential(type1=0, type2=0, potential=lj)

        self.system.addInteraction(interaction)

        interaction_energy = interaction.computeEnergy()
        # Particles does not interact across the box in x-direction
        self.assertEqual(interaction_energy, -0.30401970314257465)

    def test_moving_periodic_inner(self):
        #espressopp.PLogger.set('DomainDecompositionFree', 'TRACE')
        self.system.bc = espressopp.bc.FreeOrthorhombicBC(
            self.system.rng, (10, 10, 10), [True, True, True])
        self.system.storage = espressopp.storage.DomainDecompositionFree(self.system, self.nodeGrid, self.cellGrid)
        particle_list = [
            (1, espressopp.Real3D(4.0, 2.0, 2.0), 1.0, espressopp.Real3D(1.0, 0.0, 0.0)),
        ]
        self.system.storage.addParticles(particle_list, 'id', 'pos', 'mass', 'v')
        self.system.storage.decompose()
        integrator  = espressopp.integrator.VelocityVerlet(self.system)
        integrator.dt  = 0.1
        integrator.run(20)
        p = self.system.storage.getParticle(1)
        self.assertAlmostEqual(p.pos[0], 6.0)

    def test_moving_non_periodic_inner(self):
        #espressopp.PLogger.set('DomainDecompositionFree', 'TRACE')
        self.system.bc = espressopp.bc.FreeOrthorhombicBC(
            self.system.rng, (10, 10, 10), [False, True, True])
        self.system.storage = espressopp.storage.DomainDecompositionFree(self.system, self.nodeGrid, self.cellGrid)
        particle_list = [
            (1, espressopp.Real3D(4.0, 2.0, 2.0), 1.0, espressopp.Real3D(1.0, 0.0, 0.0)),
        ]
        self.system.storage.addParticles(particle_list, 'id', 'pos', 'mass', 'v')
        self.system.storage.decompose()
        integrator  = espressopp.integrator.VelocityVerlet(self.system)
        integrator.dt  = 0.1
        integrator.run(40)
        p = self.system.storage.getParticle(1)
        self.assertAlmostEqual(p.pos[0], 8.0)  # Cross the bounduary between nodes, nothing happen

    def test_moving_non_periodic_bounduary(self):
        #espressopp.PLogger.set('DomainDecompositionFree', 'TRACE')
        self.system.bc = espressopp.bc.FreeOrthorhombicBC(
            self.system.rng, (10, 10, 10), [False, True, True])
        self.system.storage = espressopp.storage.DomainDecompositionFree(self.system, self.nodeGrid, self.cellGrid)
        particle_list = [
            (1, espressopp.Real3D(4.0, 2.0, 2.0), 1.0, espressopp.Real3D(1.0, 0.0, 0.0)),
        ]
        self.system.storage.addParticles(particle_list, 'id', 'pos', 'mass', 'v')
        self.system.storage.decompose()

        integrator  = espressopp.integrator.VelocityVerlet(self.system)
        integrator.dt  = 0.01

        for _ in range(50):
            integrator.run(10)
        p = self.system.storage.getParticle(1)
        self.assertAlmostEqual(p.pos[0], 9.0)

    def test_moving_periodic_bounduary_x(self):
        """Test if particle can freely move across the nodes."""
        self.system.bc = espressopp.bc.FreeOrthorhombicBC(
            self.system.rng, (10, 10, 10), [True, True, True])
        self.system.storage = espressopp.storage.DomainDecompositionFree(self.system, self.nodeGrid, self.cellGrid)
        particle_list = [
            (1, espressopp.Real3D(4.0, 2.0, 2.0), 1.0, espressopp.Real3D(1.0, 0.0, 0.0)),
        ]
        self.system.storage.addParticles(particle_list, 'id', 'pos', 'mass', 'v')
        self.system.storage.decompose()

        integrator  = espressopp.integrator.VelocityVerlet(self.system)
        integrator.dt  = 0.01

        for _ in range(70):
            p = self.system.storage.getParticle(1)
            integrator.run(10)
        self.assertAlmostEqual(p.pos[0], (70*10*0.01+4.0) % 10.0)

    def test_moving_periodic_bounduary_y(self):
        """Test if particle can freely move across the nodes."""
        self.system.bc = espressopp.bc.FreeOrthorhombicBC(
            self.system.rng, (10, 10, 10), [True, True, True])
        self.system.storage = espressopp.storage.DomainDecompositionFree(self.system, self.nodeGrid, self.cellGrid)
        particle_list = [
            (1, espressopp.Real3D(4.0, 4.0, 2.0), 1.0, espressopp.Real3D(0.0, 1.0, 0.0)),
        ]
        self.system.storage.addParticles(particle_list, 'id', 'pos', 'mass', 'v')
        self.system.storage.decompose()

        integrator  = espressopp.integrator.VelocityVerlet(self.system)
        integrator.dt  = 0.01

        for _ in range(80):
            p = self.system.storage.getParticle(1)
            integrator.run(10)
        self.assertAlmostEqual(p.pos[1], (800*0.01+4.0) % 10.0)

    def test_moving_periodic_bounduary_z(self):
        """Test if particle can freely move across the nodes."""
        self.system.bc = espressopp.bc.FreeOrthorhombicBC(
            self.system.rng, (10, 10, 10), [True, True, True])
        self.system.storage = espressopp.storage.DomainDecompositionFree(self.system, self.nodeGrid, self.cellGrid)
        particle_list = [
            (1, espressopp.Real3D(4.0, 4.0, 4.0), 1.0, espressopp.Real3D(0.0, 0.0, 1.0)),
        ]
        self.system.storage.addParticles(particle_list, 'id', 'pos', 'mass', 'v')
        self.system.storage.decompose()

        integrator  = espressopp.integrator.VelocityVerlet(self.system)
        integrator.dt  = 0.01

        for _ in range(80):
            p = self.system.storage.getParticle(1)
            integrator.run(10)
        self.assertAlmostEqual(p.pos[2], (800*0.01+4.0) % 10.0)

if __name__ == '__main__':
    unittest.main()
