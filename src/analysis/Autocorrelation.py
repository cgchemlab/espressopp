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
*************************************
**espressopp.analysis.Autocorrelation**
*************************************

"""
from espressopp.esutil import cxxinit
from espressopp import pmi

from _espressopp import analysis_Autocorrelation

class AutocorrelationLocal(analysis_Autocorrelation):
    'The (local) storage of configurations.'
    def __init__(self, system):
      cxxinit(self, analysis_Autocorrelation, system)
    def gather(self, value):
      return self.cxxclass.gather(self, value)
    def clear(self):
      return self.cxxclass.clear(self)
      
    def compute(self):
      return self.cxxclass.compute(self)
    
if pmi.isController:
  class Autocorrelation(object):
    """Class for parallel analysis"""
    __metaclass__ = pmi.Proxy
    pmiproxydefs = dict(
      cls =  'espressopp.analysis.AutocorrelationLocal',
      pmicall = [ "gather", "clear", "compute" ],
      localcall = ["__getitem__", "all"],
      pmiproperty = ["size"]
    )
