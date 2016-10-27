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
            (2, 2, espressopp.Real3D(2.5, 2.0, 2.0), 2, 1),
            (3, 2, espressopp.Real3D(2.5, 2.0, 2.0), 2, 2),
            (4, 4, espressopp.Real3D(2.5, 2.0, 2.0), 2, 1)
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


class TestATRPActivator(ESPPTestCase):
    """Compare the result of previous implementation with current"""
    def test_wrong_state(self):
        atrp_activator = espressopp.integrator.ATRPActivator(self.system, 1, 100, 0.2, 0.8, 0.04, 0.01, 0.9)
        with self.assertRaises(Exception) as cm:
            atrp_activator.add_reactive_center(
                type_id=3,
                min_state=2,
                max_state=0,
                new_property=espressopp.ParticleProperties(type=3),
                delta_state=1)
        self.assertEqual(cm.exception.message, 'Min_state > max_state')

    def test_overlap_states(self):
        atrp_activator = espressopp.integrator.ATRPActivator(self.system, 1, 100, 0.2, 0.8, 0.04, 0.01, 0.9)
        atrp_activator.add_reactive_center(type_id=3, min_state=0, max_state=2,
                                           new_property=espressopp.ParticleProperties(type=3),
                                           delta_state=1)
        with self.assertRaises(Exception) as cm:
            atrp_activator.add_reactive_center(type_id=3, min_state=1, max_state=3,
                                               new_property=espressopp.ParticleProperties(type=3),
                                               delta_state=1)
        self.assertEqual(cm.exception.message, 'Min/max state overlaps')

    def test_non_reactive_centers(self):
        atrp_activator = espressopp.integrator.ATRPActivator(self.system, 1, 100, 0.2, 0.8, 0.04, 0.01, 0.9)
        atrp_activator.add_reactive_center(
            type_id=3,
            min_state=0,
            max_state=2,
            new_property=espressopp.ParticleProperties(type=3),
            delta_state=1)
        self.integrator.addExtension(atrp_activator)
        self.integrator.run(1000)
        # Non reactive centers defined
        self.assertEqual(self.system.storage.getParticle(1).state, 1)
        self.assertEqual(self.system.storage.getParticle(2).state, 1)

    def test_deactivate_particle(self):
        atrp_activator = espressopp.integrator.ATRPActivator(self.system, 1, 100, 0.2, 0.8, 0.04, 0.01, 0.9)
        # Define particles of type 3 that will be deactivate if it is in state \in [0,2)
        atrp_activator.add_reactive_center(
            type_id=2,
            min_state=0,
            max_state=2,
            new_property=espressopp.ParticleProperties(type=3),
            delta_state=1)
        self.integrator.addExtension(atrp_activator)
        self.integrator.run(1000)
        # Non reactive centers defined
        self.assertEqual(self.system.storage.getParticle(1).state, 1)
        self.assertEqual(self.system.storage.getParticle(2).state, 2)
        self.assertEqual(self.system.storage.getParticle(2).type, 3)  # Change type because of new property


if __name__ == '__main__':
    unittest.main()
