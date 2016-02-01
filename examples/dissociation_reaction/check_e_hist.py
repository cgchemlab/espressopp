#! /usr/bin/env python
#
# Copyright (c) 2016 Jakub Krajniak <jkrajniak@gmail.com>
#
# Distributed under terms of the GNU GPLv3 license.
# 

from matplotlib import pyplot as plt
import numpy
import sys

data = numpy.loadtxt(sys.argv[1], skiprows=2)

column = int(sys.argv[2])
c_data = data[:, column]

plt.hist(c_data, normed=True, bins=50)
plt.show()
