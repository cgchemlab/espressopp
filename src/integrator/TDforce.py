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


"""
*******************************
**espresso.integrator.TDforce**
*******************************

"""
from espresso.esutil import cxxinit
from espresso import pmi

from _espresso import integrator_TDforce

class TDforceLocal(integrator_TDforce):
    'The (local) Velocity Verlet Integrator.'
    def __init__(self, system, center=None):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, integrator_TDforce, system)

            # set center of TD force
            if center:
                self.cxxclass.setCenter(self, center[0], center[1], center[2])

    def addForce(self, itype, filename, cg_type, force_type=1):
            """Adds the TD force.

            Args:
                itype: The type of the interpolation.
                filename: The tabulated force file name.
                cg_type: The type of the CG particle that would be affected by the force.
                force_type: The direction of the force, if it should be x-direction then put 0
                  if the radius then put 1.
            """
            if pmi.workerIsActive():
                self.cxxclass.addForce(self, itype, filename, cg_type, force_type)

if pmi.isController :
    class TDforce(object):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls =  'espresso.integrator.TDforceLocal',
            pmiproperty = [ 'itype', 'filename'],
            pmicall = ['addForce']
            )
