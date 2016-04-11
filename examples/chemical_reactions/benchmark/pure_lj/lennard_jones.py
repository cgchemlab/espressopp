import espressopp

import cPickle
import os
import time


def define_particles(system, integrator, Npart):
    print "adding ", Npart, " particles to the system ..."
    particle_list = []
    particle_prop = ('id', 'pos', 'res_id', 'state')
    for pid in range(Npart):
        pos = system.bc.getRandomPos()
        particle_list.append((pid, pos, pid, 1))
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


def save_conf(system, Npart, filename='eq_conf.pck'):
    # Save coordinate file
    dump_particles = []
    for pid in range(Npart):
        p = system.storage.getParticle(pid)
        dump_particles.append((pid, p.type, p.pos, p.v, p.res_id, p.state))
    particle_prop = ('id', 'type', 'pos', 'v', 'res_id', 'state')
    cPickle.dump((particle_prop, dump_particles), open(filename, 'wb'))
    print('Saved configuration to {}'.format(filename))


# Settings
execfile('conf.py')

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
print "eq_conf            = ", eq_conf


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
has_eq = os.path.exists(eq_conf)
if has_eq:
    print('Reading EQ configuration from {}'.format(eq_conf))
    particle_prop, particle_list = cPickle.load(open(eq_conf, 'rb'))
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
system.addInteraction(interaction, 'lj')
system.storage.cellAdjust()

system.storage.decompose()

if not has_eq:
    # Run additional equilibration.
    for step in range(warmup_nloops):
        integrator.run(warmup_isteps)
    print('Finished equilibration')

integrator.resetTimers()
integrator.step = 0

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
system.storage.decompose()

if not has_eq:
    # Run additional equilibration.
    for step in range(warmup_nloops):
        integrator.run(warmup_isteps)
    print('Finished equilibration')
    save_conf(system, Npart, eq_conf)

integrator.resetTimers()
integrator.step = 0

system_analysis = espressopp.analysis.SystemMonitor(
    system,
    integrator,
    espressopp.analysis.SystemMonitorOutputCSV('energy.csv'))
temp_comp = espressopp.analysis.Temperature(system)
system_analysis.add_observable('T', temp_comp)
system_analysis.add_observable(
    'Ekin', espressopp.analysis.KineticEnergy(system, temp_comp))

for label, interaction in sorted(system.getAllInteractions().items()):
    system_analysis.add_observable(
        label,
        espressopp.analysis.PotentialEnergy(system, interaction))

system_analysis.add_observable(
    'fpl', espressopp.analysis.NFixedPairListEntries(system, fpl_a_a))

ext_analysis = espressopp.integrator.ExtAnalyze(system_analysis, 10)
integrator.addExtension(ext_analysis)

print('set up dump_h5md')
traj_file = espressopp.io.DumpH5MD(
    system,
    'traj.h5',
    group_name='atoms',
    static_box=False,
    author='xxx',
    email='xxx',
    store_species=True,
    store_state=True)

dump_topol = espressopp.io.DumpTopology(system, integrator, traj_file)
dump_topol.observe_tuple(fpl_a_a, 'fpl_a_a')

dump_topol.dump()
dump_topol.update()

ext_dump = espressopp.integrator.ExtAnalyze(dump_topol, 100)
integrator.addExtension(ext_dump)

print "starting benchmark..."
system_analysis.dump()
system_analysis.info()
traj_file.dump(integrator.step, integrator.step*dt)
t0 = time.time()
for step in range(equil_nloops):
    # perform equilibration_isteps integration steps
    integrator.run(equil_isteps)
    system_analysis.info()
    traj_file.dump(integrator.step, integrator.step*dt)
    dump_topol.update()
    if step % 10 == 0:
        traj_file.flush()
    # print status information
else:
    system_analysis.dump()
    system_analysis.info()
    dump_topol.update()
    traj_file.dump(equil_isteps*equil_nloops, equil_isteps*equil_nloops*dt)
    traj_file.close()
t1 = time.time()
time_file = open('benchmark_data.csv', 'a+')
time_file.write('{} {:e}\n'.format(NCPUs, t1 - t0))

save_conf(system, Npart, 'state.pck')
print "finished"
