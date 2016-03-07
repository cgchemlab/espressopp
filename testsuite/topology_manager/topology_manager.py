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
        topology_manager.rebuild()
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

    def test_merge_res_id(self):
        """Checks if two molecules will get the same res_id after bond is formed."""
        self.fpl3.addBonds([(3, 5)])
        self.topology_manager.exchange_data()
        self.integrator.run(1)
        res_id_after = {self.system.storage.getParticle(i).res_id for i in range(1, self.N+1)}
        assert res_id_after == set([1])

if __name__ == '__main__':
    unittest.main()
