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


"""Python functions to print timings from C++."""

import sys

def show(alltimers, system=None, precision=None):
    """Prints the timers data collected from all nodes.

    Args:
        alltimers: The timers.
        system: The system object.
    """
    skip_timers = ['timeRun']
    nprocs = len(alltimers)
    timers = {k: 0.0 for k, _ in alltimers[0]}
    alltimers = [{k: float(v) for k, v in ntimer if k not in skip_timers} for ntimer in alltimers]
    for ntimer in alltimers:
        for k, v in ntimer.items():
            timers[k] += v
    total_t = 0.0
    for k in timers:
        total_t += timers[k]
        timers[k] /= nprocs
    t = total_t / nprocs

    # There is a set of timers each for the interaction. The order is the same as the
    # interactions in the system object.
    sys.stdout.write('{:25} time {:3.6f}\n'.format('Run', t))

    if system is not None:
        for k, v in sorted(timers.items()):
            if k.startswith('f'):
                lbl = system.getNameOfInteraction(int(k.replace('f', '')))
            else:
                lbl = k
            sys.stdout.write('{:25} time {:3.6f} ({:3.0f} %)\n'.format(lbl, v, 100*v/t))

    sys.stdout.write('\n')
