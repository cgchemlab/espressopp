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


r"""
*****************************************
**espressopp.esutil.ParticlePairScaling**
*****************************************
"""

from espressopp import pmi
from espressopp.esutil import cxxinit
from _espressopp import esutil_ParticlePairScaling

class ParticlePairScalingLocal(esutil_ParticlePairScaling):
  def __init__(self, default_scale, incr_scale_factor, vl, integrator):
      if pmi.workerIsActive():
          cxxinit(self, esutil_ParticlePairScaling, default_scale, incr_scale_factor, vl, integrator)

if pmi.isController:
    class ParticlePairScaling:
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls = 'espressopp.esutil.ParticlePairScalingLocal',
        )
    
