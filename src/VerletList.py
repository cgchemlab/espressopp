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
***********************
**espressopp.VerletList**
***********************

"""
from espressopp import pmi
import _espressopp 
import espressopp
from espressopp.esutil import cxxinit

class VerletListLocal(_espressopp.VerletList):
    'The (local) verlet list.'

    def __init__(self, system, cutoff, exclusionlist=[]):
        'Local construction of a verlet list'
        if pmi.workerIsActive():
            if (exclusionlist == []):
                # rebuild list in constructor
                cxxinit(self, _espressopp.VerletList, system, cutoff, True)
            else:
                # do not rebuild list in constructor
                cxxinit(self, _espressopp.VerletList, system, cutoff, False)
                # add exclusions
                for pair in exclusionlist:
                    pid1, pid2 = pair
                    self.cxxclass.exclude(self, pid1, pid2)
                # now rebuild list with exclusions
                self.cxxclass.rebuild(self)
                
            
    def totalSize(self):
        'count number of pairs in VerletList, involves global reduction'
        if pmi.workerIsActive():
            return self.cxxclass.totalSize(self)
        
    def localSize(self):
        'count number of pairs in local VerletList'
        if pmi.workerIsActive():
            return self.cxxclass.localSize(self)
        
    def exclude(self, exclusionlist):
        """
        Each processor takes the broadcasted exclusion list
        and adds it to its list.
        """
        if pmi.workerIsActive():
            for pair in exclusionlist:
                pid1, pid2 = pair
                self.cxxclass.exclude(self, pid1, pid2)
            # rebuild list with exclusions
            self.cxxclass.rebuild(self)
            
    def getAllPairs(self):
        'return the pairs of the local verlet list'
        if pmi.workerIsActive():
            pairs=[]
            npairs=self.localSize()
            for i in range(npairs):
              pair=self.cxxclass.getPair(self, i+1)
              pairs.append(pair)
            return pairs 


if pmi.isController:
  class VerletList(object):
    __metaclass__ = pmi.Proxy
    pmiproxydefs = dict(
      cls = 'espressopp.VerletListLocal',
      pmiproperty = [ 'builds' ],
      pmicall = [ 'totalSize', 'exclude', 'connect', 'disconnect', 'getVerletCutoff' ],
      pmiinvoke = [ 'getAllPairs' ]
    )
