"""Semi-unit test for checking of ChemLab Framework."""

import espressopp  # pylint:disable=F0401
import unittest

try:
    import MPI  # pylint: disable=F0401
except ImportError:
    from mpi4py import MPI


class ESPPTestCase(unittest.TestCase):
    def setUp(self):
        box = (10, 10, 10)
        system = espressopp.System()
        system.kb = 1.0
        system.rng = espressopp.esutil.RNG()
        system.bc = espressopp.bc.OrthorhombicBC(system.rng, box)
        system.skin = 0.3
        self.system = system

        nodeGrid = espressopp.tools.decomp.nodeGrid(MPI.COMM_WORLD.size)
        cellGrid = espressopp.tools.decomp.cellGrid(box, nodeGrid, 2.5, system.skin)
        system.storage = espressopp.storage.DomainDecomposition(system, nodeGrid, cellGrid)

        self.integrator = espressopp.integrator.VelocityVerlet(system)
        self.integrator.dt = 0.0025

        vl = espressopp.VerletList(system, cutoff=1.0)

        self.part_prop = ('id', 'type', 'pos', 'res_id', 'state')
        particle_list = [
            (1, 1, espressopp.Real3D(2.0, 2.0, 2.0), 1, 1),
            (2, 2, espressopp.Real3D(2.5, 2.0, 2.0), 2, 1)
        ]
        system.storage.addParticles(particle_list, *self.part_prop)

        self.fpl1 = espressopp.FixedPairList(system.storage)

        topology_manager = espressopp.integrator.TopologyManager(system)
        topology_manager = topology_manager
        topology_manager.observe_tuple(self.fpl1)
        topology_manager.initialize_topology()
        self.topology_manager = topology_manager
        self.integrator.addExtension(topology_manager)

        self.ar = espressopp.integrator.ChemicalReaction(
            system, vl, system.storage, topology_manager, 1)
        self.integrator.addExtension(self.ar)


class TestBasicReaction(ESPPTestCase):

    def test_synthesis_reaction(self):
        """Runs synthesis reaction of two particles."""
        r_type_1 = espressopp.integrator.Reaction(
            type_1=1,
            type_2=2,
            delta_1=1,
            delta_2=1,
            min_state_1=1,
            max_state_1=3,
            min_state_2=1,
            max_state_2=3,
            rate=400.0,
            fpl=self.fpl1,
            cutoff=1.0)
        self.ar.add_reaction(r_type_1)

        fpl1_before = self.fpl1.getBonds()[0]
        self.integrator.run(10)
        fpl1_after = self.fpl1.getBonds()[0]

        assert fpl1_before != fpl1_after
        assert fpl1_after == [(1, 2)]

    def test_synthesis_reaction_wrong_state(self):
        """Runs synthesis reaction but with not valid conditions for min/max state."""
        r_type_1 = espressopp.integrator.Reaction(
            type_1=1,
            type_2=2,
            delta_1=1,
            delta_2=1,
            min_state_1=3,
            max_state_1=4,
            min_state_2=3,
            max_state_2=4,
            rate=400.0,
            fpl=self.fpl1,
            cutoff=1.0)
        self.ar.add_reaction(r_type_1)

        fpl1_before = self.fpl1.getBonds()[0]
        self.integrator.run(10)
        fpl1_after = self.fpl1.getBonds()[0]
        # The fpl should not change because reaction does not happen.
        assert fpl1_before == fpl1_after

    def test_synthesis_reaction_cutoff_small(self):
        """Runs synthesis reaction but with too short cutoff."""
        r_type_1 = espressopp.integrator.Reaction(
            type_1=1,
            type_2=2,
            delta_1=1,
            delta_2=1,
            min_state_1=1,
            max_state_1=4,
            min_state_2=1,
            max_state_2=4,
            rate=400.0,
            fpl=self.fpl1,
            cutoff=0.1)
        self.ar.add_reaction(r_type_1)

        fpl1_before = self.fpl1.getBonds()[0]
        self.integrator.run(10)
        fpl1_after = self.fpl1.getBonds()[0]
        # The fpl should not change because reaction does not happen.
        assert fpl1_before == fpl1_after

    def test_synthesis_reaction_rate_small(self):
        """Runs synthesis reaction but with rate = 0.0."""
        r_type_1 = espressopp.integrator.Reaction(
            type_1=1,
            type_2=2,
            delta_1=1,
            delta_2=1,
            min_state_1=1,
            max_state_1=4,
            min_state_2=1,
            max_state_2=4,
            rate=0.0,
            fpl=self.fpl1,
            cutoff=1.0)
        self.ar.add_reaction(r_type_1)

        fpl1_before = self.fpl1.getBonds()[0]
        self.integrator.run(10)
        fpl1_after = self.fpl1.getBonds()[0]
        # The fpl should not change because reaction does not happen.
        assert fpl1_before == fpl1_after


class TestCaseChangeNeighbourProperty(ESPPTestCase):
    def setUp(self):
        super(TestCaseChangeNeighbourProperty, self).setUp()
        particle_list = [
            (3, 3, espressopp.Real3D(3.0, 2.0, 2.0), 2, 1),
            (4, 3, espressopp.Real3D(3.5, 2.0, 2.0), 2, 1),
            (5, 3, espressopp.Real3D(3.5, 2.0, 2.0), 2, 1)
        ]
        self.system.storage.addParticles(particle_list, *self.part_prop)

        self.fpl1.addBonds([(2, 3), (3, 4), (4, 5)])
        self.topology_manager.exchange_data()

    def test_reaction_1(self):
        r_type_1 = espressopp.integrator.Reaction(
            type_1=1,
            type_2=2,
            delta_1=1,
            delta_2=1,
            min_state_1=1,
            max_state_1=4,
            min_state_2=1,
            max_state_2=4,
            rate=400.0,
            fpl=self.fpl1,
            cutoff=0.6)

        # Define post-process that will change a type of particle 4 from 4 to 5.
        pp_type_1 = espressopp.integrator.PostProcessChangeNeighboursProperty(
            self.topology_manager)
        pp_type_1.add_change_property(
            3, espressopp.ParticleProperties(5, 1.0, 0.0), 2)
        r_type_1.add_postprocess(pp_type_1, 'type_2')

        self.ar.add_reaction(r_type_1)
        self.integrator.run(10)

        # Check the types of particles.
        assert [self.system.storage.getParticle(x).type for x in range(1, 5)] == [1, 2, 3, 5]

    def test_reaction_2(self):
        r_type_1 = espressopp.integrator.Reaction(
            type_1=1,
            type_2=2,
            delta_1=1,
            delta_2=1,
            min_state_1=1,
            max_state_1=4,
            min_state_2=1,
            max_state_2=4,
            rate=400.0,
            fpl=self.fpl1,
            cutoff=0.6)

        # Define post-process that will change a type of particle 4 from 3 to 5.
        # and change particle 3 from 3 to 7 (different separation = 1!)
        pp_type_1 = espressopp.integrator.PostProcessChangeNeighboursProperty(
            self.topology_manager)
        pp_type_1.add_change_property(
            3, espressopp.ParticleProperties(5, 1.0, 0.0), 2)
        pp_type_1.add_change_property(
            3, espressopp.ParticleProperties(7, 1.0, 0.0), 1)
        r_type_1.add_postprocess(pp_type_1, 'type_2')

        self.ar.add_reaction(r_type_1)
        self.integrator.run(10)

        # Check the types of particles.
        assert [self.system.storage.getParticle(x).type for x in range(1, 6)] == [1, 2, 7, 5, 3]


class TestCaseChangePropertyOnState(ESPPTestCase):

    def test_reaction_1(self):
        r_type_1 = espressopp.integrator.Reaction(
            type_1=1,
            type_2=2,
            delta_1=1,
            delta_2=1,
            min_state_1=1,
            max_state_1=4,
            min_state_2=1,
            max_state_2=4,
            rate=400.0,
            fpl=self.fpl1,
            cutoff=0.6)
        pp_type_1 = espressopp.integrator.PostProcessChangePropertyOnState()

        # Change type of particle (current type = 1) to type 5.
        pp_type_1.add_change_property(
            1, espressopp.ParticleProperties(5, 1.0, 0.0), 2)
        # Change type of particle (current type = 1) to type 5.
        pp_type_1.add_change_property(
            2, espressopp.ParticleProperties(5, 1.0, 0.0), 2)
        r_type_1.add_postprocess(pp_type_1, 'both')

        self.ar.add_reaction(r_type_1)

        fpl1_before = self.fpl1.getBonds()[0]
        self.integrator.run(10)
        fpl1_after = self.fpl1.getBonds()[0]
        assert fpl1_after == [(1, 2)]

        p_types = [self.system.storage.getParticle(x).type for x in range(1, 3)]
        assert p_types == [5, 5]


class TestCaseMultipleNodes(ESPPTestCase):
    def test_node(self):
        particle_list = [
            (3, 1, espressopp.Real3D(2.5, 2.5, 4.9), 3, 1),
            (4, 1, espressopp.Real3D(2.5, 2.5, 5.1), 4, 1)
        ]
        self.system.storage.addParticles(particle_list, *self.part_prop)
        self.system.storage.decompose()

        print 'node p1', self.system.storage.mapPositionToNodeClipped(
            espressopp.Real3D(2.5, 2.5, 4.9))
        print 'node p2', self.system.storage.mapPositionToNodeClipped(
            espressopp.Real3D(2.5, 2.5, 5.1))
        r_type_1 = espressopp.integrator.Reaction(
            type_1=1,
            type_2=2,
            delta_1=1,
            delta_2=1,
            min_state_1=1,
            max_state_1=3,
            min_state_2=1,
            max_state_2=3,
            rate=400.0,
            fpl=self.fpl1,
            cutoff=1.0)
        r_type_2 = espressopp.integrator.Reaction(
            type_1=1,
            type_2=2,
            delta_1=1,
            delta_2=1,
            min_state_1=1,
            max_state_1=3,
            min_state_2=1,
            max_state_2=3,
            rate=400.0,
            fpl=self.fpl1,
            cutoff=1.0)
        self.ar.add_reaction(r_type_1)
        self.ar.add_reaction(r_type_2)
        fpl1_before = self.fpl1.getBonds()
        print fpl1_before
        self.integrator.run(2)
        fpl1_after = self.fpl1.getBonds()
        print fpl1_after


if __name__ == '__main__':
    unittest.main()
