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

.. function:: espressopp.FixedPairLambdaList.getDist(pid1, pid2)

		:param pid1: 
		:param pid2: 
		:type pid1: 
		:type pid2: 
		:rtype: 

.. function:: espressopp.FixedPairLambdaList.getPairs()

		:rtype: 

.. function:: espressopp.FixedPairLambdaList.getPairsDist()

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
            for bond in bondlist:
                pid1, pid2 = bond
                self.cxxclass.add(self, pid1, pid2)

    def getPairs(self):

        if pmi.workerIsActive():
          bonds=self.cxxclass.getPairs(self)
          return bonds 

    def getPairsDist(self):

        if pmi.workerIsActive():
          bonds=self.cxxclass.getPairsDist(self)
          return bonds 
        
    def getDist(self, pid1, pid2):
        if pmi.workerIsActive():
          return self.cxxclass.getDist(self, pid1, pid2)
        
if pmi.isController:
  class FixedPairLambdaList(object):
    __metaclass__ = pmi.Proxy
    pmiproxydefs = dict(
        cls = 'espressopp.FixedPairLambdaListLocal',
        localcall = [ "add" ],
        pmicall = [ "addPairs" ],
        pmiinvoke = ['getPairs', 'getPairsDist', 'size']
    )
    
    def getDist(self, pid1, pid2):
      pairs = pmi.invoke(self.pmiobject, 'getDist', pid1, pid2)
      for i in pairs:
        if( i != -1 ):
          return i
