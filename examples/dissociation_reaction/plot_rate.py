#! /usr/bin/env python
#
# Copyright (c) 2016 Jakub Krajniak <jkrajniak@gmail.com>
#
# Distributed under terms of the GNU GPLv3 license.
# 

import h5py
import argparse
import numpy as np

from matplotlib import pyplot as plt

parser = argparse.ArgumentParser()
parser.add_argument('h5')

args = parser.parse_args()

h5 = h5py.File(args.h5, 'r')

bonds = np.array(h5['/connectivity/fpl/value'])

bonds_time = [x[x != -1].shape[0] / 2 for x in bonds]

plt.plot(bonds_time)
plt.show()
