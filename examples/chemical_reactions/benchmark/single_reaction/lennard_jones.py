"""
This is an example for an MD simulation of a simple Lennard-Jones fluid
with ESPResSo++. We will start with particles at random positions within
the simulation box interacting via a shifted Lennard-Jones type potential
with an interaction cutoff at 2.5.
Newtons equations of motion are integrated with a Velocity-Verlet integrator.
The canonical (NVT) ensemble is realized by using a Langevin thermostat.
In order to prevent explosion due to strongly overlapping volumes of
random particles the system needs to be warmed up first.
Warm-up is accomplished by using a repelling-only LJ interaction
(cutoff=1.12246, shift=0.25) with a force capping at radius 0.6
and initial small LJ epsilon value of 0.1.
During warmup epsilon is gradually increased to its final value 1.0.
After warm-up the system is equilibrated using the full uncapped  LJ Potential.

If a system still explodes during warmup or equilibration, warmup time
could be increased by increasing warmup_nloops and the capradius could
be set to another value. Depending on the system (number of particles, density, ...)
it could also be necessary to vary sigma during warmup.

The simulation consists of the following steps:

  1. specification of the main simulation parameters
  2. setup of the system, random number generator and parallelisation
  3. setup of the integrator and simulation ensemble
  4. adding the particles
  5. setting up interaction potential for the warmup
  6. running the warmup loop
  7. setting up interaction potential for the equilibration
  8. running the equilibration loop
  9. writing configuration to a file
"""

import espressopp

import cPickle
import logging
import os
import time


def define_particles(system, integrator, Npart):
    print "adding ", Npart, " particles to the system ..."
    particle_list = []
    particle_prop = ('id', 'type', 'pos', 'res_id', 'state')
    for pid in range(Npart):
        pos = system.bc.getRandomPos()
        particle_list.append((pid, 0, pos, pid, 2))
    system.storage.addParticles(particle_list, *particle_prop)
    system.storage.decompose()

    verletlist = espressopp.VerletList(system, warmup_cutoff)
    LJpot = espressopp.interaction.LennardJonesCapped(
        epsilon=epsilon_start,
        sigma=sigma,
        cutoff=warmup_cutoff,
        caprad=capradius,
        shift='auto')
    interaction = espressopp.interaction.VerletListLennardJonesCapped(verletlist)
    interaction.setPotential(type1=0, type2=0, potential=LJpot)

    system.addInteraction(interaction)
    print "starting warm-up ..."
    espressopp.tools.info(system, integrator)
    for step in range(warmup_nloops):
        integrator.run(warmup_isteps)
        LJpot.epsilon = epsilon_delta
        interaction.setPotential(type1=0, type2=0, potential=LJpot)
        espressopp.tools.info(system, integrator)
    print "warmup finished"
    system.removeInteraction(0)
    verletlist.disconnect()

    # Save coordinate file
    dump_particles = []
    for pid in range(Npart):
        p = system.storage.getParticle(pid)
        dump_particles.append((pid, p.type, p.pos, p.v, p.res_id, p.state))
    particle_prop = ('id', 'type', 'pos', 'v', 'res_id', 'state')
    cPickle.dump((particle_prop, dump_particles), open('eq_conf.pck', 'wb'))


# Settings
Npart              = 3272
rho                = 0.8442
L                  = pow(Npart/rho, 1.0/3.0)
box                = (L, L, L)
r_cutoff           = 2.5
skin               = 1.5
temperature        = 1.0
dt                 = 0.005
epsilon            = 1.0
sigma              = 1.0

warmup_cutoff      = pow(2.0, 1.0/6.0)
warmup_nloops      = 100
warmup_isteps      = 200
total_warmup_steps = warmup_nloops * warmup_isteps
epsilon_start      = 0.1
epsilon_end        = 1.0
epsilon_delta      = (epsilon_end - epsilon_start) / warmup_nloops
capradius          = 0.6
equil_nloops       = 100
equil_isteps       = 100

print espressopp.Version().info()
print "Npart              = ", Npart
print "rho                = ", rho
print "L                  = ", L
print "box                = ", box
print "r_cutoff           = ", r_cutoff
print "skin               = ", skin
print "temperature        = ", temperature
print "dt                 = ", dt
print "epsilon            = ", epsilon
print "sigma              = ", sigma
print "warmup_cutoff      = ", warmup_cutoff
print "warmup_nloops      = ", warmup_nloops
print "warmup_isteps      = ", warmup_isteps
print "total_warmup_steps = ", total_warmup_steps
print "epsilon_start      = ", epsilon_start
print "epsilon_end        = ", epsilon_end
print "epsilon_delta      = ", epsilon_delta
print "capradius          = ", capradius
print "equil_nloops       = ", equil_nloops
print "equil_isteps       = ", equil_isteps


system             = espressopp.System()
system.rng         = espressopp.esutil.RNG(12345)
system.bc          = espressopp.bc.OrthorhombicBC(system.rng, box)
system.skin        = skin
NCPUs              = espressopp.MPI.COMM_WORLD.size
nodeGrid           = espressopp.tools.decomp.nodeGrid(NCPUs)
cellGrid           = espressopp.tools.decomp.cellGrid(box, nodeGrid, warmup_cutoff, skin)
system.storage     = espressopp.storage.DomainDecomposition(system, nodeGrid, cellGrid)

print "NCPUs              = ", NCPUs
print "nodeGrid           = ", nodeGrid
print "cellGrid           = ", cellGrid

integrator     = espressopp.integrator.VelocityVerlet(system)
integrator.dt  = dt
if (temperature != None):
  # create e Langevin thermostat
  thermostat             = espressopp.integrator.LangevinThermostat(system)
  # set Langevin friction constant
  thermostat.gamma       = 1.0
  # set temperature
  thermostat.temperature = temperature
  # tell the integrator to use this thermostat
  integrator.addExtension(thermostat)

# Read particles or equilibrate it.
has_eq = os.path.exists('eq_conf.pck')
if has_eq:
    print('Reading EQ configuration from eq_conf.pck')
    particle_prop, particle_list = cPickle.load(open('eq_conf.pck', 'rb'))
    system.storage.addParticles(particle_list, *particle_prop)
else:
    print('EQ configuration not found, running equilibration process for N={}'.format(Npart))
    define_particles(system, integrator, Npart)

# Set up interactions
dynamic_ex_list = espressopp.DynamicExcludeList(integrator)
verletlist = espressopp.VerletList(system, r_cutoff, exclusionlist=dynamic_ex_list)
interaction = espressopp.interaction.VerletListLennardJones(verletlist)
potential = interaction.setPotential(type1=0, type2=0,
                                     potential=espressopp.interaction.LennardJones(
                                         epsilon=epsilon, sigma=sigma, cutoff=r_cutoff, shift=0.0))
system.addInteraction(interaction)
system.storage.cellAdjust()

system.storage.decompose()

if not has_eq:
    # Run additional equilibration.
    for step in range(warmup_nloops):
        integrator.run(warmup_isteps)
    print('Finished equilibration')

integrator.resetTimers()
integrator.step = 0

print('Defining a chemical reaction')
for i in range(Npart):
    system.storage.modifyParticle(i, 'state', 2)

fpl_a_a = espressopp.FixedPairList(system.storage)
potHarmonic = espressopp.interaction.Harmonic(
    K=30.0,
    r0=0.97,
    cutoff=1.5)
interHarmonic = espressopp.interaction.FixedPairListHarmonic(
    system,
    fpl_a_a,
    potHarmonic)
fpl_a_a.addBonds([])
system.addInteraction(interHarmonic, 'harmonic')
dynamic_ex_list.observe_tuple(fpl_a_a)

ar = espressopp.integrator.ChemicalReaction(
    system,
    verletlist,
    system.storage,
    10)

r_type_1 = espressopp.integrator.Reaction(
    type_1=0,
    type_2=0,
    delta_1=-1,
    delta_2=-1,
    min_state_1=1,
    max_state_1=3,
    min_state_2=1,
    max_state_2=3,
    rate=10.0,
    fpl=fpl_a_a,
    cutoff=1.5)

r_type_1.intramolecular=True
ar.add_reaction(r_type_1)
integrator.addExtension(ar)

#topology_manager = espressopp.integrator.TopologyManager(system)
#topology_manager.rebuild()
#topology_manager.observe_tuple(fpl_a_a)
#topology_manager.initialize_topology()
#integrator.addExtension(topology_manager)

#logging.getLogger('Reaction').setLevel(logging.DEBUG)
#logging.getLogger('ChemicalReaction').setLevel(logging.DEBUG)

print "starting benchmark..."
t0 = time.time()
for step in range(equil_nloops):
    print step
    # perform equilibration_isteps integration steps
    integrator.run(equil_isteps)
    # print status information
t1 = time.time()
time_file = open('benchmark_data.csv', 'a+')
time_file.write('{:e}\n'.format(t1 - t0))
print "finished"
print "has: ", fpl_a_a.totalSize()
