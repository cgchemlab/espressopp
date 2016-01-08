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
********************************
**espressopp.FixedPairLambdaList**
********************************


.. function:: espressopp.FixedPairLambdaList(storage, initial_lambda)

		:param storage: 
		:type storage:
		:param initial_lambda:
		:type initial_lambda: float

.. function:: espressopp.FixedPairLambdaList.add(pid1, pid2)

		:param pid1: 
		:param pid2: 
		:type pid1: 
		:type pid2: 
		:rtype: 

.. function:: espressopp.FixedPairLambdaList.addPairs(bondlist)

		:param bondlist: 
		:type bondlist: 
		:rtype: 

.. function:: espressopp.FixedPairLambdaList.getLambda(pid1, pid2)

		:param pid1: 
		:param pid2: 
		:type pid1: 
		:type pid2: 
		:rtype: 

.. function:: espressopp.FixedPairLambdaList.getPairs()

		:rtype: 

.. function:: espressopp.FixedPairLambdaList.getPairsLambda()

		:rtype: 

.. function:: espressopp.FixedPairLambdaList.size()

		:rtype: 
"""
from espressopp import pmi
import _espressopp 
import espressopp
from espressopp.esutil import cxxinit

class FixedPairLambdaListLocal(_espressopp.FixedPairLambdaList):
    def __init__(self, storage, initial_lambda=1.0):
        if pmi.workerIsActive():
            cxxinit(self, _espressopp.FixedPairLambdaList, storage, initial_lambda)

    def add(self, pid1, pid2):
        if pmi.workerIsActive():
            return self.cxxclass.add(self, pid1, pid2)

    def size(self):
        if pmi.workerIsActive():
            return self.cxxclass.size(self)

    def addPairs(self, bondlist):
        """
        Each processor takes the broadcasted bondlist and
        adds those pairs whose first particle is owned by
        this processor.
        """
        if pmi.workerIsActive():
            for pid1, pid2 in bondlist:
                self.cxxclass.add(self, pid1, pid2)

    def getPairs(self):
        if pmi.workerIsActive():
          return self.cxxclass.getPairs(self)

    def getPairsDist(self):
        if pmi.workerIsActive():
          return self.cxxclass.getPairsDist(self)
        
    def getLambda(self, pid1, pid2):
        if pmi.workerIsActive():
          return self.cxxclass.getLambda(self, pid1, pid2)
        
if pmi.isController:
  class FixedPairLambdaList(object):
      __metaclass__ = pmi.Proxy
      pmiproxydefs = dict(
          cls = 'espressopp.FixedPairLambdaListLocal',
          localcall = [ "add" ],
          pmicall = [ "addPairs" ],
          pmiinvoke = ['getPairs', 'getPairsLambda', 'size', 'getLambda']
      )