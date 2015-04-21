#  Copyright (C) 2015
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


"""
************************************
**espressopp.integrator.FixDistances**
************************************

"""
from espressopp.esutil import cxxinit
from espressopp import pmi
from espressopp.integrator.Extension import *
from _espressopp import integrator_FixDistances

class FixDistancesLocal(ExtensionLocal, integrator_FixDistances):
    'The (local) Fix Positions part.'
    def __init__(self, system, cs_list=None, anchor_type=None, target_type=None):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            if anchor_type is not None and target_type is not None:
                cxxinit(self, integrator_FixDistances, system, anchor_type, target_type)
            else:
                cxxinit(self, integrator_FixDistances, system)
            if cs_list:
                self.addConstraints(cs_list)

    def addConstraints(self, cs_list):
        for anchor_id, target_id, dist in cs_list:
            self.cxxclass.add_triplet(self, anchor_id, target_id, dist)


if pmi.isController :
    class FixDistances(Extension):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls =  'espressopp.integrator.FixDistancesLocal',
            pmicall = ['addConstraints', 'add_postprocess'],
            pmiproperty = ['size']
            )
