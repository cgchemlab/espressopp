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
**********************************
**espressopp.FixedPairListLambda**
**********************************


.. function:: espressopp.FixedPairListLambda(storage, initial_lambda)

		:param storage: 
		:type storage:
		:param initial_lambda:
		:type initial_lambda: float

.. function:: espressopp.FixedPairListLambda.add(pid1, pid2)

		:param pid1: 
		:param pid2: 
		:type pid1: 
		:type pid2: 
		:rtype: 

.. function:: espressopp.FixedPairListLambda.addPairs(bondlist)

		:param bondlist: 
		:type bondlist: 
		:rtype: 

.. function:: espressopp.FixedPairListLambda.getLambda(pid1, pid2)

		:param pid1: 
		:param pid2: 
		:type pid1: 
		:type pid2: 
		:rtype: 

.. function:: espressopp.FixedPairListLambda.getBonds()

		:rtype: 

.. function:: espressopp.FixedPairListLambda.getPairsLambda()

		:rtype: 

.. function:: espressopp.FixedPairListLambda.size()

		:rtype: 
"""
from espressopp import pmi
import _espressopp 
from espressopp.esutil import cxxinit


class FixedPairListLambdaLocal(_espressopp.FixedPairListLambda):
    def __init__(self, storage, initial_lambda=1.0):
        if pmi.workerIsActive():
            cxxinit(self, _espressopp.FixedPairListLambda, storage, initial_lambda)

    def add(self, pid1, pid2):
        if pmi.workerIsActive():
            return self.cxxclass.add(self, pid1, pid2)

    def size(self):
        if pmi.workerIsActive():
            return self.cxxclass.size(self)

    def addBonds(self, bondlist):
        if pmi.workerIsActive():
            for bond in bondlist:
                pid1, pid2 = bond
                self.cxxclass.add(self, pid1, pid2)

    def getPairs(self):
        if pmi.workerIsActive():
          return self.cxxclass.getPairs(self)
        
    def getLambda(self, pid1, pid2):
        if pmi.workerIsActive():
          return self.cxxclass.getLambda(self, pid1, pid2)


if pmi.isController:
  class FixedPairListLambda(object):
      __metaclass__ = pmi.Proxy
      pmiproxydefs = dict(
          cls = 'espressopp.FixedPairListLambdaLocal',
          localcall = [ "add" ],
          pmicall = [ "addBonds", "setLambda", "setAllLambda"],
          pmiinvoke = ['getBonds', 'getPairsLambda', 'size', 'getLambda']
      )