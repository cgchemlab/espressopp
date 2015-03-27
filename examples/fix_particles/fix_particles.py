import espressopp

# create default Lennard Jones (WCA) system with 0 particles and cubic box (L=10)
system, integrator = espressopp.standard_system.LennardJones(0, (10*1.12, 10*1.12, 10*1.12))

C_FIXED = 1
C_FREE  = 0
# fix x,y and z coord axis
fixMask = espressopp.Int3D(C_FREE, C_FIXED, C_FREE)

# create a particel group that will contain the fixed particles
fixedWall  = espressopp.ParticleGroup(system.storage)

# add a particle wall
pid = 1
for k in range(10):
  for l in range(10):
    system.storage.addParticle(pid, espressopp.Real3D(k*1.12, 5, l*1.12))
    fixedWall.add(pid)
    pid += 1

# add also one free particle
system.storage.addParticle(0, espressopp.Real3D(5.8,9,5.5))
system.storage.modifyParticle(0, 'v', espressopp.Real3D(0, -0.1, 0))

# don't forget do decompose !
system.storage.decompose()

# create FixPositions Extension and add it to the integrator
fixpositions = espressopp.integrator.FixPositions(system, fixedWall, fixMask)
integrator.addExtension(fixpositions)

# run the simulation
sock = espressopp.tools.vmd.connect(system)
for i in range(10000):
  integrator.run(100)
  espressopp.tools.vmd.imd_positions(system, sock)
