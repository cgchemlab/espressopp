#  Copyright (c) 2015
#      Jakub Krajniak (jkrajniak at gmail.com)
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


from espressopp.esutil import cxxinit
from espressopp import pmi

from espressopp.integrator.Extension import *
from _espressopp import integrator_DynamicResolution
from _espressopp import integrator_BasicDynamicResolutionType


class BasicDynamicResolutionLocal(ExtensionLocal, integrator_BasicDynamicResolutionType):
    """The (local) BasicDynamicResolutionType

    Args:
        system: The system object
        type_rate: The dictionary with type id and rate.
    """
    def __init__(self, system, type_rate):
        'Local constructor'
        if pmi.workerIsActive():
            cxxinit(self, integrator_BasicDynamicResolution)
            for type_id, rate in type_rate.iteritems():
                self.cxxclass.set_type_rate(type_id, rate)

    def set_type_rate(type_id, rate):
        if pmi.workerIsActive():
            self.cxxclass.set_type_rate(type_id, rate)


class DynamicResolutionLocal(ExtensionLocal, integrator_DynamicResolution):
    'The (local) DynamicResolution'

    def __init__(self, _system, _fixedtuplelist, _rate):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, integrator_DynamicResolution, _system,  _fixedtuplelist, _rate)


if pmi.isController:
    class BasicDynamicResolution(Extension):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = {
            'cls': 'espressopp.integrator.BasicDynamicResolutionLocal',
            'pmicall': ['set_type_rate']
        }

    class DynamicResolution(Extension):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls = 'espressopp.integrator.DynamicResolutionLocal',
            pmiproperty = [ 'resolution', 'rate', 'active' ],
            pmicall = ['SetPosVel']
            )