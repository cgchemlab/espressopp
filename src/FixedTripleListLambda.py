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
************************************
**espressopp.FixedTripleListLambda**
************************************


.. function:: espressopp.FixedTripleListLambda(storage)

		:param storage:
		:type storage:

.. function:: espressopp.FixedTripleListLambda.add(pid1, pid2, pid3)

		:param pid1:
		:param pid2:
		:param pid3:
		:type pid1:
		:type pid2:
		:type pid3:
		:rtype:

.. function:: espressopp.FixedTripleListLambda.addTriples(triplelist)

		:param triplelist:
		:type triplelist:
		:rtype:

.. function:: espressopp.FixedTripleListLambda.getTriples()

		:rtype:

.. function:: espressopp.FixedTripleListLambda.size()

		:rtype:
"""
from espressopp import pmi
import _espressopp
import espressopp
from espressopp.esutil import cxxinit


class FixedTripleListLambdaLocal(_espressopp.FixedTripleListLambda):
    def __init__(self, storage, initial_lambda=1.0):
        if pmi.workerIsActive():
            cxxinit(self, _espressopp.FixedTripleListLambda, storage, initial_lambda)

    def add(self, pid1, pid2, pid3):
        if pmi.workerIsActive():
            return self.cxxclass.add(self, pid1, pid2, pid3)

    def addTriples(self, triplelist):
        if pmi.workerIsActive():
            for triple in triplelist:
                pid1, pid2, pid3 = triple
                self.cxxclass.add(self, pid1, pid2, pid3)

    def size(self):
        if pmi.workerIsActive():
            return self.cxxclass.size(self)

    def getTriples(self):
        if pmi.workerIsActive():
            triples = self.cxxclass.getTriples(self)
            return triples

    def getLambda(self, pid1, pid2, pid3):
        if pmi.workerIsActive():
            return self.cxxclass.getLambda(self, pid1, pid2, pid3)


if pmi.isController:
    class FixedTripleListLambda(object):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls = 'espressopp.FixedTripleListLambdaLocal',
            localcall = [ 'add' ],
            pmicall = [ 'addTriples', 'totalSize', 'setLambda', 'setLambdaAll'],
            pmiinvoke = ['getTriples', 'size', 'getLambda', 'getPairsLambda']
        )