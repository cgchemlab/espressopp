"""Semi-unit test for checking of Topology Manager."""

import espressopp  # pylint:disable=F0401
import unittest

try:
    import MPI  # pylint: disable=F0401
except ImportError:
    from mpi4py import MPI


class TestTopologyManager(unittest.TestCase):
    def setUp(self):
        # Initialize the espressopp system
        box = (10, 10, 10)
        system = espressopp.System()
        self.system = system
        system.kb = 1.0
        system.rng = espressopp.esutil.RNG()
        system.bc = espressopp.bc.OrthorhombicBC(system.rng, box)
        system.skin = 0.3

        nodeGrid = espressopp.tools.decomp.nodeGrid(MPI.COMM_WORLD.size)
        cellGrid = espressopp.tools.decomp.cellGrid(box, nodeGrid, 2.5, system.skin)
        system.storage = espressopp.storage.DomainDecomposition(system, nodeGrid, cellGrid)

        # Adding ten nodes
        self.N = 8
        particle_list = []
        for pid in range(1, self.N+1):
            pos = system.bc.getRandomPos()
            particle_list.append((pid, pos, 1 if pid < 5 else 2))
        system.storage.addParticles(particle_list, 'id', 'pos', 'res_id')
        system.storage.decompose()
        self.integrator = espressopp.integrator.VelocityVerlet(system)
        self.integrator.dt = 0.0025

        self.fpl1 = espressopp.FixedPairList(system.storage)
        self.fpl1.addBonds([(1, 2), (2, 3), (2, 4)])
        self.fpl2 = espressopp.FixedPairList(system.storage)
        self.fpl2.addBonds([(5, 6), (6, 7), (6, 8)])
        self.fpl3 = espressopp.FixedPairList(system.storage)

        self.ftl = espressopp.FixedTripleList(system.storage)
        self.ftl.addTriples([(1, 2, 3), (1, 2, 4), (5, 6, 7), (5, 6, 8)])

        self.fql = espressopp.FixedQuadrupleList(system.storage)
        self.fql2 = espressopp.FixedQuadrupleList(system.storage)

        topology_manager = espressopp.integrator.TopologyManager(system)
        self.topology_manager = topology_manager
        topology_manager.observe_tuple(self.fpl1)
        topology_manager.observe_tuple(self.fpl2)
        topology_manager.observe_tuple(self.fpl3)
        topology_manager.register_triplet(self.ftl, 0)
        topology_manager.register_quadruplet(self.fql, 0)
        topology_manager.register_quadruplet(self.fql2, 0, 0, 0, 1)
        topology_manager.register_tuple(self.fpl3, 0, 0)
        topology_manager.initialize_topology()
        self.integrator.addExtension(topology_manager)

    def test_topology(self):
        """Tests if adjacent list is built correctly."""
        nb_list = dict(self.topology_manager.get_neighbour_lists()[0])
        # Check nb
        assert nb_list[1] == [2]
        assert set(nb_list[2]) == set([1, 4, 3])
        assert set(nb_list[3]) == set([2])
        assert set(nb_list[4]) == set([2])
        assert set(nb_list[5]) == set([6])
        assert set(nb_list[6]) == set([5, 7, 8])
        assert set(nb_list[7]) == set([6])
        assert set(nb_list[8]) == set([6])

    def test_add_new_bond(self):
        """Test if topology was updated after adding new bond."""
        self.fpl3.addBonds([(3, 5)])
        self.topology_manager.exchange_data()
        nb_list = dict(self.topology_manager.get_neighbour_lists()[0])
        assert nb_list[1] == [2]
        assert set(nb_list[2]) == set([1, 4, 3])
        assert set(nb_list[3]) == set([2, 5])
        assert set(nb_list[4]) == set([2])
        assert set(nb_list[5]) == set([3, 6])
        assert set(nb_list[6]) == set([5, 7, 8])
        assert set(nb_list[7]) == set([6])
        assert set(nb_list[8]) == set([6])

    def test_triplets(self):
        """Test content of triplets"""
        ftl_before = set(self.ftl.getTriples()[0])
        self.fpl3.addBonds([(3, 5)])
        self.topology_manager.exchange_data()
        ftl_after = set(self.ftl.getTriples()[0])
        # Checks if triplet list is extended correctly
        assert ftl_before != ftl_after
        assert not ((ftl_after - ftl_before) - set([(3, 5, 6), (2, 3, 5), (5, 3, 2), (6, 5, 3)]))

    def test_quadruplets(self):
        """Test the content of quadruplets."""
        self.fpl3.addBonds([(3, 5)])
        self.topology_manager.exchange_data()
        # Checks quadruplets.
        new_quadruples = set(self.fql.getQuadruples()[0])
        assert not (new_quadruples - set([(5, 3, 2, 1), (5, 3, 2, 4), (3, 5, 6, 7),
                                          (3, 5, 6, 8), (2, 3, 5, 6)]))
        # Check if second quadruple is empty as types does not match.
        assert not self.fql2.getQuadruples()[0]

    def test_residues_bonds(self):
        """Test if residues are updated correctly."""
        self.fpl3.addBonds([(3, 5)])
        self.topology_manager.exchange_data()
        assert self.topology_manager.is_residue_connected(1,2)

    def test_molecules_split(self):
        self.assertEqual(len(self.topology_manager.get_molecule_ids()[0]), 2)
        self.fpl3.addBonds([(3, 5)])
        self.topology_manager.exchange_data()
        self.assertEqual(len(self.topology_manager.get_molecule_ids()[0]), 1)
        self.fpl3.remove(3, 5, False)
        self.topology_manager.exchange_data()
        self.assertEqual(len(self.topology_manager.get_molecule_ids()[0]), 2)


class TestTopologyCycle(unittest.TestCase):
    def setUp(self):
        # Initialize the espressopp system
        box = (10, 10, 10)
        system = espressopp.System()
        self.system = system
        system.kb = 1.0
        system.rng = espressopp.esutil.RNG()
        system.bc = espressopp.bc.OrthorhombicBC(system.rng, box)
        system.skin = 0.3

        nodeGrid = espressopp.tools.decomp.nodeGrid(MPI.COMM_WORLD.size)
        cellGrid = espressopp.tools.decomp.cellGrid(box, nodeGrid, 2.5, system.skin)
        system.storage = espressopp.storage.DomainDecomposition(system, nodeGrid, cellGrid)

        # Adding ten nodes
        self.N = 8
        particle_list = []
        for pid in range(0, self.N):
            pos = system.bc.getRandomPos()
            particle_list.append((pid+1, pos, pid/2+1))
        system.storage.addParticles(particle_list, 'id', 'pos', 'res_id')
        system.storage.decompose()
        self.integrator = espressopp.integrator.VelocityVerlet(system)
        self.integrator.dt = 0.0025

        self.fpl1 = espressopp.FixedPairList(system.storage)
        self.fpl1.addBonds([(i, i+1) for i in range(1, self.N+1, 2)])

        self.fpl3 = espressopp.FixedPairList(system.storage)

        topology_manager = espressopp.integrator.TopologyManager(system)
        self.topology_manager = topology_manager
        topology_manager.observe_tuple(self.fpl1)
        topology_manager.observe_tuple(self.fpl3)
        topology_manager.register_tuple(self.fpl3, 0, 0)
        topology_manager.initialize_topology()
        self.integrator.addExtension(topology_manager)

    def test_connect_two_molecules(self):
        self.assertEqual(self.topology_manager.get_molecule_ids()[0], [1, 2, 3, 4])
        self.fpl3.addBonds([(2, 3)])
        self.topology_manager.exchange_data()
        self.assertEqual(self.topology_manager.get_molecule_ids()[0], [1, 3, 4])
        self.assertEqual(self.topology_manager.get_molecule(1)[0], [1, 2, 3, 4])
        self.assertTrue(self.topology_manager.is_particle_connected(2, 3))
        self.assertTrue(self.topology_manager.is_particle_connected(3, 4))
        self.assertFalse(self.topology_manager.is_particle_connected(4, 5))

    def test_connect_into_ring(self):
        self.assertEqual(self.topology_manager.get_molecule_ids()[0], [1, 2, 3, 4])
        self.fpl3.addBonds([(2, 3)])
        self.topology_manager.exchange_data()
        self.assertEqual(self.topology_manager.get_molecule_ids()[0], [1, 3, 4])
        self.fpl3.addBonds([(4, 5)])
        self.topology_manager.exchange_data()
        self.assertEqual(self.topology_manager.get_molecule_ids()[0], [1, 4])
        self.assertEqual(self.topology_manager.get_molecule(1)[0], [1, 2, 3, 4, 5, 6])
        self.fpl3.addBonds([(6, 7)])
        self.topology_manager.exchange_data()
        self.assertEqual(self.topology_manager.get_molecule_ids()[0], [1])
        self.assertEqual(self.topology_manager.get_molecule(1)[0], [1, 2, 3, 4, 5, 6, 7, 8])
        self.topology_manager.exchange_data()
        self.fpl3.addBonds([(1, 8)])
        self.topology_manager.exchange_data()
        self.assertEqual(self.topology_manager.get_molecule_ids()[0], [1])
        self.assertEqual(self.topology_manager.get_molecule(1)[0], [1, 2, 3, 4, 5, 6, 7, 8])
        self.topology_manager.exchange_data()
        self.assertTrue(self.topology_manager.is_particle_connected(1, 8))

    def test_break_ring(self):
        self.assertEqual(self.topology_manager.get_molecule_ids()[0], [1, 2, 3, 4])
        self.fpl3.addBonds([(2, 3), (4, 5), (6, 7), (1, 8)])
        self.topology_manager.exchange_data()
        self.assertEqual(self.topology_manager.get_molecule_ids()[0], [1])
        self.assertEqual(self.topology_manager.get_molecule(1)[0], [1, 2, 3, 4, 5, 6, 7, 8])
        # Break one bond, still one molecule
        self.assertTrue(self.topology_manager.is_particle_connected(4, 5))
        self.fpl3.remove(4, 5)
        self.topology_manager.exchange_data()
        self.assertEqual(self.topology_manager.get_molecule_ids()[0], [1])
        self.assertEqual(self.topology_manager.get_molecule(1)[0], [1, 2, 3, 4, 5, 6, 7, 8])
        self.assertFalse(self.topology_manager.is_particle_connected(4, 5))
        # Break second bond, now two molecules
        self.assertTrue(self.topology_manager.is_particle_connected(6, 7))
        self.fpl3.remove(6, 7)
        self.topology_manager.exchange_data()
        self.assertEqual(self.topology_manager.get_molecule_ids()[0], [1, 5])
        self.fpl3.addBonds([(6, 7)])
        self.topology_manager.exchange_data()
        self.assertEqual(self.topology_manager.get_molecule_ids()[0], [5])


class TestCheckPropertyAtDistance(unittest.TestCase):
    def setUp(self):
        # Initialize the espressopp system
        box = (10, 10, 10)
        system = espressopp.System()
        self.system = system
        system.kb = 1.0
        system.rng = espressopp.esutil.RNG()
        system.bc = espressopp.bc.OrthorhombicBC(system.rng, box)
        system.skin = 0.3

        nodeGrid = espressopp.tools.decomp.nodeGrid(MPI.COMM_WORLD.size)
        cellGrid = espressopp.tools.decomp.cellGrid(box, nodeGrid, 2.5, system.skin)
        system.storage = espressopp.storage.DomainDecomposition(system, nodeGrid, cellGrid)

        # Adding ten nodes
        self.N = 6
        particle_prop = ['id', 'pos', 'type', 'res_id', 'mass', 'state']
        particle_list = [
            [1, None, 1, 1, 1.0, 2],
            [2, None, 1, 1, 1.0, 2],
            [3, None, 3, 1, 1.0, 3],
            [4, None, 2, 1, 1.0, 9],
            [5, None, 4, 1, 1.0, 8],
            [6, None, 3, 1, 1.0, 7]
        ]
        for pid in range(0, self.N):
            pos = system.bc.getRandomPos()
            particle_list[pid][1] = pos
        system.storage.addParticles(particle_list, *particle_prop)
        system.storage.decompose()
        self.integrator = espressopp.integrator.VelocityVerlet(system)
        self.integrator.dt = 0.0025

        bonds = [(1, 2), (1, 5), (1, 6), (2, 3), (2, 4)]

        self.fpl1 = espressopp.FixedPairList(system.storage)
        self.fpl1.addBonds(bonds)

        self.fpl3 = espressopp.FixedPairList(system.storage)

        topology_manager = espressopp.integrator.TopologyManager(system)
        self.topology_manager = topology_manager
        topology_manager.observe_tuple(self.fpl1)
        topology_manager.observe_tuple(self.fpl3)
        topology_manager.register_tuple(self.fpl3, 0, 0)
        topology_manager.initialize_topology()
        self.integrator.addExtension(topology_manager)

    def test_has_property_neighbour(self):
        pp = espressopp.integrator.TopologyParticleProperties(3)
        pp.set_min_max_state(7, 8)
        self.assertTrue(self.topology_manager.has_neighbour_particle_property(3, pp, 3))

if __name__ == '__main__':
    unittest.main()
