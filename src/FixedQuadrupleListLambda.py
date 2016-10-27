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
***************************************
**espressopp.FixedQuadrupleListLambda**
***************************************


.. function:: espressopp.FixedQuadrupleListLambda(storage)

		:param storage: 
		:type storage: 

.. function:: espressopp.FixedQuadrupleListLambda.add(pid1, pid2, pid3, pid4)

		:param pid1: 
		:param pid2: 
		:param pid3: 
		:param pid4: 
		:type pid1: 
		:type pid2: 
		:type pid3: 
		:type pid4: 
		:rtype: 

.. function:: espressopp.FixedQuadrupleListLambda.addQuadruples(quadruplelist)

		:param quadruplelist: 
		:type quadruplelist: 
		:rtype: 

.. function:: espressopp.FixedQuadrupleListLambda.getQuadruples()

		:rtype: 

.. function:: espressopp.FixedQuadrupleListLambda.size()

		:rtype: 
"""
import _espressopp
from espressopp import pmi
from espressopp.esutil import cxxinit


class FixedQuadrupleListLambdaLocal(_espressopp.FixedQuadrupleListLambda):
    def __init__(self, storage, init_lambda=1.0):
        if pmi.workerIsActive():
            cxxinit(self, _espressopp.FixedQuadrupleListLambda, storage, init_lambda)

    def add(self, pid1, pid2, pid3, pid4):
        if pmi.workerIsActive():
            return self.cxxclass.add(self, pid1, pid2, pid3, pid4)

    def size(self):
        if pmi.workerIsActive():
            return self.cxxclass.size(self)

    def addQuadruples(self, quadruplelist):
        if pmi.workerIsActive():
            for quadruple in quadruplelist:
                pid1, pid2, pid3, pid4 = quadruple
                self.cxxclass.add(self, pid1, pid2, pid3, pid4)

    def getQuadruples(self):
        if pmi.workerIsActive():
            quadruple = self.cxxclass.getQuadruples(self)
            return quadruple

    def getLambda(self, pid1, pid2, pid3, pid4):
        if pmi.workerIsActive():
            return self.cxxclass.getLambda(self, pid1, pid2, pid3, pid4)

if pmi.isController:
    class FixedQuadrupleListLambda(object):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls='espressopp.FixedQuadrupleListLambdaLocal',
            localcall=["add"],
            pmicall=["addQuadruples", "totalSize", 'setLambda', 'setAllLambda'],
            pmiinvoke=["getQuadruples", "size", 'getLambda', 'getPairsLambda']
        )
