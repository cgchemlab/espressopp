#! /usr/bin/env python
#
# Copyright (c) 2015 Jakub Krajniak <jkrajniak@gmail.com>
#
# Distributed under terms of the GNU GPLv3 license.
#

import espressopp
try:
    import MPI  # pylint: disable=F0401
except ImportError:
    from mpi4py import MPI

system = espressopp.System()
system.rng = espressopp.esutil.RNG()
system.bc = espressopp.bc.OrthorhombicBC(system.rng, (10.0, 10.0, 10.0))
system.skin = 0.3

nodeGrid = espressopp.tools.decomp.nodeGrid(MPI.COMM_WORLD.size)
cellGrid = espressopp.tools.decomp.cellGrid((10.0, 10.0, 10.0), nodeGrid, 1.5, 0.3)
print(nodeGrid, cellGrid)
system.storage = espressopp.storage.DomainDecomposition(system, nodeGrid, cellGrid)

# Generate 100 particles
for i in range(100):
    pos = system.bc.getRandomPos()
    system.storage.addParticle(i+1, pos)

system.storage.decompose()

fpl = espressopp.FixedPairList(system.storage)
ftl = espressopp.FixedTripleList(system.storage)
fql = espressopp.FixedQuadrupleList(system.storage)

topology_manager = espressopp.integrator.TopologyManager(system)
topology_manager.observe(fpl)
topology_manager.observe_triple(ftl, 0, 0, 0)
topology_manager.observe_quadruple(fql, 0, 0, 0, 0)

import logging
logging.getLogger("TopologyManager").setLevel(logging.DEBUG)

fpl.addBonds([(1, 2)])
fpl.addBonds([(1, 3)])
fpl.addBonds([(2, 3)])
fpl.addBonds([(3, 2)])
fpl.addBonds([(1, 4), (3, 4)])
topology_manager.print_topology()

print('-------TRIPLETS-----------')
print ftl.getTriples()

print('-------QUADRUPLETS--------')
print fql.getQuadruples()
