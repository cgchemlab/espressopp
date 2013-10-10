# DEMONSTRATION OF THE LATTICE-BOLTZMANN SIMULATION
#
# This is a script for a simple test employing LB model. It incorporates an
# initial system in three dimensions specified by the box size. The liquid has
# a uniform density (equal to 1. at every lattice site) and is subject to a 
# harmonic force in z-direction as a function of the lateral position x, i.e.
# f_z = A * sin (2 /pi x / L_x). The amplitude A is chosen to be 0.0005 to pro-
# vide low Mach number, i.e. f_z|max = 0.0005 that is smaller than the speed of
# sound (1./3.).
import espresso
import cProfile, pstats
from espresso import Int3D
from espresso import Real3D

# create default Lennard Jones (WCA) system with 0 particles and cubic box (L=40)
system, integrator = espresso.standard_system.LennardJones(0, box=(40, 40, 40), temperature=1.0)

# define a LB grid
lb = espresso.integrator.LatticeBoltzmann(system, Ni=Int3D(40, 40, 40))

# declare gammas responsible for viscosities (if they differ from 0)
lb.gamma_b = 0.5
lb.gamma_s = 0.5

# specify desired temperature (set the fluctuations if any)
lb.lbTemp = 0.0
#lb.lbTemp = 0.00001

# ask to set the external force to zero. The program then goes to the function initExtForce
# and sets the force to a harmonic one.
lb.extForce = Real3D(0.,0.,0.)

# add extension to the integrator
integrator.addExtension(lb)

# add some profiling statistics for the run
cProfile.run("integrator.run(200)",'profiler_stats')
p = pstats.Stats('profiler_stats')
p.strip_dirs().sort_stats("time").print_stats(10)
