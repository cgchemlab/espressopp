#  Copyright (C) 2016
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
**********************************
**espressopp.analysis.Resolution**
**********************************

"""
from espressopp.esutil import cxxinit
from espressopp import pmi

from espressopp.analysis.Observable import *
from _espressopp import analysis_Resolution
from _espressopp import analysis_ResolutionFixedPairList

class ResolutionLocal(ObservableLocal, analysis_Resolution):
    def __init__(self, system):
        if pmi.workerIsActive():
            cxxinit(self, analysis_Resolution, system)

class ResolutionFixedPairListLocal(ObservableLocal, analysis_ResolutionFixedPairList):
    def __init__(self, system, fpl):
        if pmi.workerIsActive():
            cxxinit(self, analysis_ResolutionFixedPairList, system, fpl)

if pmi.isController :
    class Resolution(Observable):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls =  'espressopp.analysis.ResolutionLocal',
            pmiproperty = ['value']
        )

    class ResolutionFixedPairList(Observable):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls =  'espressopp.analysis.ResolutionFixedPairListLocal',
            pmiproperty = ['value']
        )