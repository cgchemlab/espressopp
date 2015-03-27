#!/usr/bin/env python
# -*- coding: utf-8 -*-

# relevant imports
import sys
import time
import espressopp
import mpi4py.MPI as MPI

import Tetracryst # Preparation of tetrahedral crystal and constuctions of bonds in tetrahedral liquid

from espressopp import Real3D, Int3D
from espressopp.tools import decomp
from espressopp.tools import timers

# integration steps, cutoff, skin, AdResS specifications
steps = 10000
timestep = 0.0005
intervals = 1000

rc = 4.5 # cutoff coarse-grained potential
rca = 1.122462048309373 # cutoff atomistic potential (cutoff (2^(1/6)), WCA)
skin = 0.4

# Parameters for the thermostat
#gamma = 2.0
#temp = 1.0

# Parameters for size of AdResS dimensions
ex_size = 5.0
hy_size = 5.0

# read equilibrated configuration file
pid, type, x, y, z, vx, vy, vz, Lx, Ly, Lz = espressopp.tools.readxyz("equilibrated_conf.xyz")

# Table for coarse-grained potential
tabCG = "table_potential.dat"

# number of CG particles
num_particlesCG = len(x)/4

# number of AT particles
num_particles = len(x)

# set up the system
sys.stdout.write('Setting up simulation ...\n')
density = num_particles / (Lx * Ly * Lz)
size = (Lx, Ly, Lz)

system = espressopp.System()
system.rng = espressopp.esutil.RNG()
system.bc = espressopp.bc.OrthorhombicBC(system.rng, size)
system.skin = skin

comm = MPI.COMM_WORLD
nodeGrid = decomp.nodeGrid(comm.size)
cellGrid = decomp.cellGrid(size, nodeGrid, rc, skin)

# (H-)AdResS domain decomposition
system.storage = espressopp.storage.DomainDecompositionAdress(system, nodeGrid, cellGrid)


# prepare AT particles
allParticlesAT = []
allParticles = []
tuples = []
for pidAT in range(num_particles):
    allParticlesAT.append([pidAT, # add here these particles just temporarily
                         Real3D(x[pidAT], y[pidAT], z[pidAT]), # position
                         Real3D(vx[pidAT], vy[pidAT], vz[pidAT]), # velocity
                         Real3D(0, 0, 0), # force
                         1, 1.0, 1]) # type, mass, is AT particle

# create CG particles
for pidCG in range(num_particlesCG):
    # we put CG molecule in first atom, later CG molecules will be positioned in the center
    cmp = espressopp.tools.AdressSetCG(4, pidCG, allParticlesAT)
    # Preparation of tuples (tuples define, which atoms belong to which CG molecules)
    tmptuple = [pidCG+num_particles]
    for pidAT2 in range(4):
        pid = pidCG*4+pidAT2
        tmptuple.append(pid)
    
    # append CG particles    
    allParticles.append([pidCG+num_particles, # CG particle has to be added first!
                         Real3D(cmp[0], cmp[1], cmp[2]), # pos
                         Real3D(0, 0, 0), # vel
                         Real3D(0, 0, 0), # force
                         0, 4.0, 0]) # type, mass, is not AT particle
    # append AT particles
    for pidAT in range(4): 
        pid = pidCG*4+pidAT
        allParticles.append([pid, # now the AT particles can be added
                            (allParticlesAT[pid])[1], # pos
                            (allParticlesAT[pid])[2], # vel
                            (allParticlesAT[pid])[3], # force
                            (allParticlesAT[pid])[4], # type
                            (allParticlesAT[pid])[5], # mass
                            (allParticlesAT[pid])[6]]) # is AT particle 
    # append tuple to tuplelist    
    tuples.append(tmptuple)
    

# add particles to system
system.storage.addParticles(allParticles, "id", "pos", "v", "f", "type", "mass", "adrat")

# create FixedTupleList object
ftpl = espressopp.FixedTupleListAdress(system.storage)

# and add the tuples
ftpl.addTuples(tuples)
system.storage.setFixedTuplesAdress(ftpl)

# add bonds between AT particles
fpl = espressopp.FixedPairListAdress(system.storage, ftpl)
bonds = Tetracryst.makebonds(len(x))
fpl.addBonds(bonds)

# decompose after adding tuples and bonds
print "Added tuples and bonds, decomposing now ..." 
system.storage.decompose()
print "done decomposing"

# AdResS Verlet list
vl = espressopp.VerletListAdress(system, cutoff=rc, adrcut=rc,
                                dEx=ex_size, dHy=hy_size,
                                adrCenter=[Lx/2, Ly/2, Lz/2])

# non-bonded potentials
# LJ Capped WCA between AT and tabulated potential between CG particles
interNB = espressopp.interaction.VerletListHadressLennardJones(vl, ftpl) # Here we need specific (H-)AdResS interaction type
potWCA  = espressopp.interaction.LennardJones(epsilon=1.0, sigma=1.0, shift='auto', cutoff=rca)
potCG = espressopp.interaction.Tabulated(itype=3, filename=tabCG, cutoff=rc) # CG
interNB.setPotentialAT(type1=1, type2=1, potential=potWCA) # AT
interNB.setPotentialCG(type1=0, type2=0, potential=potCG) # CG
system.addInteraction(interNB)

# bonded potentials
# Quartic potential between AT particles
potQuartic = espressopp.interaction.Quartic(K=75.0, r0=1.0)
interQuartic = espressopp.interaction.FixedPairListQuartic(system, fpl, potQuartic)
system.addInteraction(interQuartic)

# VelocityVerlet integrator
integrator = espressopp.integrator.VelocityVerlet(system)
integrator.dt = timestep

# add AdResS extension
adress = espressopp.integrator.Adress(system, vl, ftpl)
integrator.addExtension(adress)

# add Langevin thermostat extension
#langevin = espressopp.integrator.LangevinThermostat(system)
#langevin.gamma = gamma
#langevin.temperature = temp
#langevin.adress = True # enable AdResS!
#integrator.addExtension(langevin)

# distribute atoms and CG molecules according to AdResS domain decomposition, place CG molecules in the center of mass 
espressopp.tools.AdressDecomp(system, integrator)

# system information
print ''
print 'AdResS Center =', [Lx/2, Ly/2, Lz/2]
print 'number of AT particles =', num_particles
print 'number of CG particles =', num_particlesCG
print 'density = %.4f' % (density)
print 'rc =', rc
print 'dt =', integrator.dt
print 'skin =', system.skin
print 'steps =', steps
print 'NodeGrid = %s' % (nodeGrid,)
print 'CellGrid = %s' % (cellGrid,)
print ''

# analysis
temperature = espressopp.analysis.Temperature(system)

fmt = '%5d %8.4f %12.3f %12.3f %12.3f %12.3f\n'
T = temperature.compute()
Ek = 0.5 * T * (3 * num_particles)
Ep = interNB.computeEnergy()
Eb = interQuartic.computeEnergy()
sys.stdout.write(' step    Temp       etotal      enonbonded    ebonded     ekinetic\n')
sys.stdout.write(fmt % (0, T, Ek + Ep + Eb, Ep, Eb, Ek))

# Timer, Steps
start_time = time.clock()
nsteps = steps / intervals

# write the start configuration to trajectory pdb-file
dump_conf_gro = espressopp.io.DumpGRO(system, integrator, filename='trajCG.gro')
dump_conf_gro_adr = espressopp.io.DumpGROAdress(system, ftpl, integrator, filename='trajAT.gro')

# integration and on the fly analysis
for s in range(1, intervals + 1):
  integrator.run(nsteps)
  step = nsteps * s
  T = temperature.compute()
  Ek = 0.5 * T * (3 * num_particles)
  Ep = interNB.computeEnergy()
  Eb = interQuartic.computeEnergy()
  sys.stdout.write(fmt % (step, T, Ek + Ep + Eb, Ep, Eb, Ek))
  dump_conf_gro.dump()
  dump_conf_gro_adr.dump()




# simulation information
end_time = time.clock()
timers.show(integrator.getTimers(), precision=3)
sys.stdout.write('Total # of neighbors = %d\n' % vl.totalSize())
sys.stdout.write('Ave neighs/atom = %.1f\n' % (vl.totalSize() / float(num_particles)))
sys.stdout.write('Neighbor list builds = %d\n' % vl.builds)
sys.stdout.write('Integration steps = %d\n' % integrator.step)
sys.stdout.write('CPU time = %.1f\n' % (end_time - start_time))

