#  Copyright (C) 2012,2013
#      Max Planck Institute for Polymer Research
#  Copyright (C) 2008,2009,2010,2011
#      Max-Planck-Institute for Polymer Research & Fraunhofer SCAI
#  
#  This file is part of ESPResSo++.
#  
#  ESPResSo++ is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#  
#  ESPResSo++ is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#  
#  You should have received a copy of the GNU General Public License
#  along with this program.  If not, see <http://www.gnu.org/licenses/>. 


import datetime
import sys
import espressopp

def info(system, integrator, per_atom=False, valid_types=None):
  NPart_comp  = espressopp.analysis.NPart(system)
  T_comp = espressopp.analysis.Temperature(system)
  if valid_types is not None:
      for type_id in valid_types:
          T_comp.add_type(type_id)
          NPart_comp.add_type(type_id)
  NPart = NPart_comp.compute()
  T      = T_comp.compute()
  P      = espressopp.analysis.Pressure(system).compute()
  Pij    = espressopp.analysis.PressureTensor(system).compute()
  step   = integrator.step
  Ek     = (3.0/2.0) * NPart * T
  Epot   = []
  Etotal = 0.0
  if per_atom:
    tot    = '%5d %10.4f %10.6f %10.6f %12.8f' % (step, T, P, Pij[3], Ek/NPart)
  else:
    tot    = '%5d %10.4f %10.6f %10.6f %12.3f' % (step, T, P, Pij[3], Ek)      
  tt     = ''
  for k in xrange(system.getNumberOfInteractions()):
    e       = system.getInteraction(k).computeEnergy()
    Etotal += e
    if per_atom:
      tot    += ' %12.8f' % (e/NPart)
      tt     += '     e%i/N    ' % k
    else:
      tot    += ' %12.3f' % e
      tt     += '      e%i     ' % k

  if per_atom:
    tot += ' %12.8f' % (Etotal/NPart + Ek/NPart)
    tt  += '   etotal/N  '
  else:
    tot += ' %12.3f' % (Etotal + Ek)
    tt  += '    etotal   '
  tot += ' %12.8f\n' % system.bc.boxL[0]
  tt  += '    boxL     \n'
  if step == 0:
    if per_atom:
      sys.stdout.write(' step      T          P        Pxy         ekin/N  N    ' + tt)
    else:
      sys.stdout.write(' step      T          P        Pxy          ekin   N    ' + tt)        
  sys.stdout.write(tot)

def final_info(system, integrator, vl, start_time, end_time):
  NPart  = espressopp.analysis.NPart(system).compute()
  espressopp.tools.timers.show(integrator.getTimers(), system)
  sys.stdout.write('Total # of neighbors = %d\n' % vl.totalSize())
  sys.stdout.write('Ave neighs/atom = %.1f\n' % (vl.totalSize() / float(NPart)))
  sys.stdout.write('Neighbor list builds = %d\n' % vl.builds)
  sys.stdout.write('Integration steps = %d\n' % integrator.step)
  sys.stdout.write('CPUs = %i CPU time per CPU = %.5f\n' % (espressopp.MPI.COMM_WORLD.size,
                                                            (end_time - start_time)/float(espressopp.MPI.COMM_WORLD.size)))

