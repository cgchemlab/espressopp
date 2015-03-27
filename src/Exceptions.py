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
**espressopp.Exceptions**
***********************

"""
import sys, traceback

class Error(Exception):
    """Raised to show unrecoverable espressopp errors.
    """
    def __init__(self, msg):
        try:
            raise Exception
        except:
            file, lineno, module, line = traceback.extract_stack()[0]
            self.msg = 'ERROR while executing ' + str(file) + ' line ' + str(lineno) + ': ' + str(line) + '\n-> ' + msg
    def __str__(self) :
        return self.msg
    def __repr__(self) :
        return str(self)

class ParticleDoesNotExistHere(Exception):
    """ Raised to indicate, that a certain Particle does not exist on a CPU
    """
    def __init__(self, msg):
        try:
            raise Exception
        except:
            self.msg = msg
    def __str__(self) :
        return self.msg
    def __repr__(self) :
        return str(self)

class UnknownParticleProperty(Exception):
    """ Raised to indicate, that a certain Particle property does not exists 
    """
    def __init__(self, msg):
        try:
            raise Exception
        except:
            self.msg = msg
    def __str__(self) :
        return self.msg
    def __repr__(self) :
        return str(self)

class MissingFixedPairList(Exception):
    """ Raised to indicate, that a FixedPairList object is missing
    """
    def __init__(self, msg):
        try:
            raise Exception
        except:
            self.msg = msg
    def __str__(self) :
        return self.msg
    def __repr__(self) :
        return str(self)
