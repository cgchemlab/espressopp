"""Semi-unit test for checking of ChemLab Framework."""

import espressopp  # pylint:disable=F0401
import logging
import unittest

try:
    import MPI  # pylint: disable=F0401
except ImportError:
    from mpi4py import MPI


class ESPPTestCase(unittest.TestCase):
    def setUp(self):
        self.system, self.integrator = self.create_system()
        self.vl = espressopp.VerletList(self.system, cutoff=2.5)
        self.part_prop = ('id', 'type', 'pos', 'res_id', 'state')
        particle_list = [
            (1, 1, espressopp.Real3D(2.0, 2.0, 2.0), 1, 1),
            (2, 2, espressopp.Real3D(2.5, 2.0, 2.0), 2, 1)
        ]
        self.system.storage.addParticles(particle_list, *self.part_prop)
        self.fpl1 = espressopp.FixedPairList(self.system.storage)

        topology_manager = espressopp.integrator.TopologyManager(self.system)
        topology_manager.observe_tuple(self.fpl1)
        topology_manager.initialize_topology()
        self.topology_manager = topology_manager
        self.integrator.addExtension(topology_manager)

        self.ar = espressopp.integrator.ChemicalReaction(
            self.system, self.vl, self.system.storage, topology_manager, 1)
        self.integrator.addExtension(self.ar)
        super(ESPPTestCase, self).setUp()

    def create_system(self):
        box = (10, 10, 10)
        system = espressopp.System()
        system.kb = 1.0
        system.rng = espressopp.esutil.RNG(12345)
        system.bc = espressopp.bc.OrthorhombicBC(system.rng, box)
        system.skin = 0.3

        nodeGrid = espressopp.tools.decomp.nodeGrid(MPI.COMM_WORLD.size)
        cellGrid = espressopp.tools.decomp.cellGrid(box, nodeGrid, 2.5, system.skin)
        system.storage = espressopp.storage.DomainDecomposition(system, nodeGrid, cellGrid)

        integrator = espressopp.integrator.VelocityVerlet(system)
        integrator.dt = 0.0025

        return system, integrator


class TestCmpAssociatieReaction(ESPPTestCase):
    """Compare the result of previous implementation with current"""
    def test_synthesis_reaction(self):
        particles = [
            (3, 3, espressopp.Real3D(4.0, 2.0, 2.0), 1, 1),
            (4, 3, espressopp.Real3D(4.0, 3.0, 2.0), 2, 0),
            (5, 3, espressopp.Real3D(4.0, 6.0, 2.0), 3, 0),
        ]
        self.system.storage.addParticles(particles, *self.part_prop)
        self.system.storage.decompose()
        fpl = espressopp.FixedPairList(self.system.storage)
        ar = espressopp.integrator.AssociationReaction(self.system, self.vl, fpl, self.system.storage)
        ar.rate = 1000.0
        ar.interval = 1
        ar.cutoff = 1.2
        ar.typeA = 3
        ar.typeB = 3
        ar.deltaA = 1
        ar.deltaB = 1
        ar.stateAMin = 0
        self.integrator.addExtension(ar)
        self.assertEqual(fpl.getAllBonds(), [])
        self.integrator.run(2000)
        self.assertEqual(fpl.getAllBonds(), [(3, 4)])

    def test_random_box_ar(self):
        system1, integrator1 = self.create_system()
        system2, integrator2 = self.create_system()

        for pid in range(3, 1000, 1):
            pos = system1.bc.getRandomPos()
            system1.storage.addParticle(pid, pos)
            system1.storage.modifyParticle(pid, 'type', 3)
            system1.storage.modifyParticle(pid, 'state', 0)
            system1.storage.modifyParticle(pid, 'res_id', pid)
            system2.storage.addParticle(pid, pos)
            system2.storage.modifyParticle(pid, 'type', 3)
            system2.storage.modifyParticle(pid, 'state', 0)
            system2.storage.modifyParticle(pid, 'res_id', pid)

        system1.storage.decompose()
        system2.storage.decompose()

        fpl_old = espressopp.FixedPairList(system1.storage)
        vl_old = espressopp.VerletList(system1, cutoff=2.5)

        ar = espressopp.integrator.AssociationReaction(system1, vl_old, fpl_old, system1.storage)
        ar.rate = 1000.0
        ar.interval = 1
        ar.cutoff = 1.2
        ar.typeA = 3
        ar.typeB = 3
        ar.deltaA = 1
        ar.deltaB = 1
        ar.stateAMin = 0
        integrator1.addExtension(ar)

        # Second system
        vl = espressopp.VerletList(system2, cutoff=2.5)
        fpl1 = espressopp.FixedPairList(system2.storage)
        topology_manager = espressopp.integrator.TopologyManager(system2)
        ar2 = espressopp.integrator.ChemicalReaction(
            system2, vl, system2.storage, topology_manager, 1)
        integrator2.addExtension(ar2)
        r_type_1 = espressopp.integrator.Reaction(
            type_1=3,
            type_2=3,
            delta_1=1,
            delta_2=1,
            min_state_1=0,
            max_state_1=10**8,
            min_state_2=0,
            max_state_2=1,
            rate=1000.0,
            cutoff=1.2,
            fpl=fpl1)
        ar2.interval = 1
        ar2.add_reaction(r_type_1)
        ar2.nearest_mode = True
        ar2.pair_distances_filename = 'pairs_distances.txt'
        #self.assertEqual(fpl_old.getAllBonds(), [])
        #self.assertEqual(fpl1.getAllBonds(), [])
        integrator1.run(1)
        integrator2.run(1)
        ar2.pair_distances_filename = ''
        integrator2.run(1)
        ar2.pair_distances_filename = 'abc.txt'
        integrator2.run(1)
        s_old = set([tuple(sorted(x)) for x in fpl_old.getAllBonds()])
        s_new = set([tuple(sorted(x)) for x in fpl1.getAllBonds()])
        #self.assertEqual(len(s_old), len(s_new))
        #self.assertItemsEqual(fpl_old.getAllBonds(), fpl1.getAllBonds())


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

        self.assertFalse(self.topology_manager.is_particle_connected(1, 2))
        self.integrator.run(10)
        fpl1_after = self.fpl1.getBonds()[0]
        self.assertItemsEqual(fpl1_after, [(1, 2)])
        self.assertTrue(self.topology_manager.is_particle_connected(1, 2))

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
            (4, 4, espressopp.Real3D(3.5, 2.0, 2.0), 2, 1),
            (5, 3, espressopp.Real3D(3.5, 2.0, 2.0), 1, 1),
            (6, 4, espressopp.Real3D(3.5, 2.0, 2.0), 1, 1)
        ]
        self.system.storage.addParticles(particle_list, *self.part_prop)

        self.fpl1.addBonds([(2, 3), (3, 4), (1, 5), (5, 6)])
        self.topology_manager.initialize_topology()

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
            3, espressopp.integrator.TopologyParticleProperties(5, 1.0, 0.0, state=7), 1)
        r_type_1.add_postprocess(pp_type_1, 'both')

        self.ar.add_reaction(r_type_1)
        self.integrator.run(10)
        print [self.system.storage.getParticle(x).state for x in range(1, 7)]

        # Check the types of particles.
        self.assertItemsEqual([self.system.storage.getParticle(x).type for x in range(1, 7)], [1, 2, 5, 4, 5, 4])

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
            3, espressopp.integrator.TopologyParticleProperties(5, 1.0, 0.0), 2)
        pp_type_1.add_change_property(
            3, espressopp.integrator.TopologyParticleProperties(7, 1.0, 0.0), 1)
        r_type_1.add_postprocess(pp_type_1, 'type_2')

        self.ar.add_reaction(r_type_1)
        self.integrator.run(10)

        # Check the types of particles.
        self.assertItemsEqual([self.system.storage.getParticle(x).type for x in range(1, 6)], [1, 2, 7, 5, 3])


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
            1, espressopp.integrator.TopologyParticleProperties(5, 1.0, 0.0), 2)
        # Change type of particle (current type = 1) to type 5.
        pp_type_1.add_change_property(
            2, espressopp.integrator.TopologyParticleProperties(5, 1.0, 0.0), 2)
        r_type_1.add_postprocess(pp_type_1, 'both')

        self.ar.add_reaction(r_type_1)

        fpl1_before = self.fpl1.getBonds()[0]
        self.integrator.run(10)
        fpl1_after = self.fpl1.getBonds()[0]
        assert fpl1_after == [(1, 2)]

        p_types = [self.system.storage.getParticle(x).type for x in range(1, 3)]
        assert p_types == [5, 5]


class TestCyclization(ESPPTestCase):
    def test_cyclization(self):
        """Check if it is possible to make an cycle"""
        particle_list = [
            (3, 2, espressopp.Real3D(2.0, 3.0, 2.0), 3, 1),
            (4, 1, espressopp.Real3D(2.5, 3.0, 2.0), 4, 1),
        ]
        self.system.storage.addParticles(particle_list, *self.part_prop)
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
            cutoff=1.1)
        self.ar.add_reaction(r_type_1)
        fpl1_before = self.fpl1.getBonds()
        self.assertEquals(fpl1_before, [[]])
        self.integrator.run(2)
        fpl1_after = self.fpl1.getBonds()
        # Full cycle.
        self.assertItemsEqual(fpl1_after[0], [(2, 4), (3, 4), (1, 2), (1, 3)])

    def test_non_cyclization(self):
        """Check if it is possible to make an cycle"""
        particle_list = [
            (3, 2, espressopp.Real3D(2.0, 3.0, 2.0), 1, 1),
            (4, 1, espressopp.Real3D(2.5, 3.0, 2.0), 1, 1),
        ]
        self.system.storage.addParticles(particle_list, *self.part_prop)
        self.system.storage.modifyParticle(1, 'res_id', 2)
        self.system.storage.modifyParticle(2, 'res_id', 2)

        self.fpl1.addBonds([(1, 2), (3, 4), (1, 3)])
        self.topology_manager.initialize_topology()

        self.system.storage.decompose()
        self.integrator.run(0)

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
            cutoff=1.1)
        r_type_1.intraresidual = False
        self.ar.add_reaction(r_type_1)
        fpl1_before = self.fpl1.getBonds()
        self.integrator.run(2)
        fpl1_after = self.fpl1.getBonds()
        # Full cycle.
        self.assertEquals(set(fpl1_after[0]), {(1, 2), (1, 3), (3, 4)})


class TestCaseRemoveNeighbourBond(ESPPTestCase):
    def setUp(self):
        super(TestCaseRemoveNeighbourBond, self).setUp()
        particle_list = [
            (3, 3, espressopp.Real3D(3.0, 2.0, 2.0), 2, 1),
            (4, 4, espressopp.Real3D(3.5, 2.0, 2.0), 2, 1),
            (5, 4, espressopp.Real3D(3.5, 2.0, 2.0), 2, 1),
            (6, 3, espressopp.Real3D(3.5, 2.0, 2.0), 2, 1),
            (7, 4, espressopp.Real3D(3.5, 2.0, 2.0), 2, 1),
        ]
        self.system.storage.addParticles(particle_list, *self.part_prop)
        self.system.storage.decompose()

        self.fpl23 = espressopp.FixedPairList(self.system.storage)
        self.fpl24 = espressopp.FixedPairList(self.system.storage)
        self.fpl34 = espressopp.FixedPairList(self.system.storage)

        self.ftl234 = espressopp.FixedTripleList(self.system.storage)
        self.ftl234.addTriples([(2, 3, 4)])

        self.fpl23.addBonds([(2, 3)])
        self.fpl24.addBonds([(2, 7)])
        self.fpl34.addBonds([(3, 4), (3, 5), (4, 6)])

        self.topology_manager.register_tuple(self.fpl24, 2, 4)
        self.topology_manager.register_tuple(self.fpl34, 3, 4)
        self.topology_manager.register_tuple(self.fpl23, 2, 3)
        self.topology_manager.register_triplet(self.ftl234, 2, 3, 4)

        self.topology_manager.initialize_topology()


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
        pp_type_1 = espressopp.integrator.PostProcessRemoveNeighbourBond(self.topology_manager)
        pp_type_1.add_bond_to_remove(2, 1, 2, 3)
        pp_type_1.add_bond_to_remove(2, 3, 3, 4)
        r_type_1.add_postprocess(pp_type_1, 'type_2')

        self.ar.add_reaction(r_type_1)
        self.integrator.run(10)

        fpl_after = self.fpl1.getAllBonds()
        self.assertEqual(fpl_after, [(1, 2)])

        # The bond 2-3 should be removed but rest should be okey
        fpl_after = self.fpl23.getAllBonds()
        self.assertEqual(fpl_after, [])
        self.assertEqual(self.fpl24.getAllBonds(), [(2, 7)])
        self.assertEqual(self.fpl34.getAllBonds(), [(3, 4), (3, 5)])
        self.assertEqual(self.ftl234.getTriples(), [[]])  # Clean up also triplets because the bond 2-3 is removed.


class TestCaseExchangeReaction(ESPPTestCase):
    def setUp(self):
        super(TestCaseExchangeReaction, self).setUp()
        import logging
        # 0-1-2-3-2-1-2-3-2-1-0
        particle_list = [
            (3, 0, espressopp.Real3D(3.0, 2.0, 2.0), 2, 1),
            (4, 1, espressopp.Real3D(3.2, 2.0, 2.0), 2, 0),
            (5, 2, espressopp.Real3D(3.4, 2.0, 2.0), 2, 2),

            (6, 3, espressopp.Real3D(3.6, 2.0, 2.0), 3, 3),

            (7, 2, espressopp.Real3D(3.8, 2.0, 2.0), 4, 2),
            (8, 1, espressopp.Real3D(4.0, 2.0, 2.0), 4, 0),
            (9, 2, espressopp.Real3D(4.2, 2.0, 2.0), 4, 2),

            (10, 3, espressopp.Real3D(4.4, 2.0, 2.0), 5, 3),

            (11, 2, espressopp.Real3D(4.6, 2.0, 2.0), 6, 2),
            (12, 1, espressopp.Real3D(4.8, 2.0, 2.0), 6, 0),
            (13, 0, espressopp.Real3D(5.0, 2.0, 2.0), 6, 2),

            (14, 4, espressopp.Real3D(4.2, 1.8, 2.0), 7, 4),
        ]
        self.system.storage.addParticles(particle_list, *self.part_prop)
        self.system.storage.decompose()
        # ('id', 'type', 'pos', 'res_id', 'state')

        self.fpl_static = espressopp.FixedPairList(self.system.storage)
        self.ftl_static = espressopp.FixedTripleList(self.system.storage)
        self.ftl123 = espressopp.FixedTripleList(self.system.storage)
        self.ftl232 = espressopp.FixedTripleList(self.system.storage)
        self.fpl_chem = espressopp.FixedPairList(self.system.storage)

        self.dynamic_exclude = espressopp.DynamicExcludeList(self.integrator)
        self.dynamic_exclude.observe_tuple(self.fpl_static)
        self.dynamic_exclude.observe_tuple(self.fpl_chem)
        self.dynamic_exclude.observe_triple(self.ftl_static)
        self.dynamic_exclude.observe_triple(self.ftl123)
        self.dynamic_exclude.observe_triple(self.ftl232)
        self.dynamic_exclude.observe_tuple(self.fpl_chem)

        self.fpl_static.addBonds([(3, 4), (4, 5), (7, 8), (8, 9), (11, 12), (12, 13)])
        self.topology_manager.observe_tuple(self.fpl_static)

        self.ftl_static.addTriples([(3, 4, 5), (7, 8, 9), (11, 12, 13)])

        self.ftl123.addTriples([(4, 5, 6), (6, 7, 8), (8, 9, 10), (10, 11, 12)])
        self.topology_manager.register_triplet(self.ftl123, 1, 2, 3)

        self.ftl232.addTriples([(5, 6, 7), (9, 10, 11)])
        self.topology_manager.register_triplet(self.ftl232, 2, 3, 2)

        self.fpl_chem.addBonds([(5, 6), (6, 7), (9, 10), (10, 11)])
        self.topology_manager.register_tuple(self.fpl_chem, 2, 3)
        self.topology_manager.initialize_topology()
        self.dynamic_exclude.update()

    def test_remove_bond(self):
        r_type_1 = espressopp.integrator.Reaction(
            type_1=2,
            type_2=4,
            delta_1=1,
            delta_2=1,
            min_state_1=2,
            max_state_1=4,
            min_state_2=4,
            max_state_2=5,
            rate=400.0,
            fpl=self.fpl_chem,
            cutoff=0.2)
        r_type_1.is_virtual = True
        pp_type_1 = espressopp.integrator.PostProcessRemoveNeighbourBond(self.topology_manager)
        pp_type_1.add_bond_to_remove(2, 1, 2, 3)
        r_type_1.add_postprocess(pp_type_1, 'type_1')

        fd = espressopp.integrator.FixDistances(self.system, None, 2, 4)
        self.integrator.addExtension(fd)

        self.assertEqual(fd.totalSize(), 0)

        pp_join_fd = espressopp.integrator.PostProcessJoinParticles(fd, 0.5)
        r_type_1.add_postprocess(pp_join_fd, 'type_1')

        self.ar.add_reaction(r_type_1)

        static_angles = self.ftl_static.getTriples()

        self.assertEqual(self.dynamic_exclude.size, 38)

        self.integrator.run(2)
        self.assertEqual(fd.totalSize(), 1)

        self.assertItemsEqual(self.fpl_chem.getAllBonds(), [(5, 6), (6, 7), (10, 11)])
        self.assertItemsEqual(self.ftl_static.getTriples(), static_angles)
        self.assertItemsEqual(self.ftl232.getAllTriples(), [(5, 6, 7)])
        self.assertItemsEqual(self.ftl123.getAllTriples(), [(4, 5, 6), (6, 7, 8), (10, 11, 12)])

        self.assertEqual(self.dynamic_exclude.size, 32)


if __name__ == '__main__':
    unittest.main()
