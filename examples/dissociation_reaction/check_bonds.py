#! /usr/bin/env python
#
# Copyright (c) 2016 Jakub Krajniak <jkrajniak@gmail.com>
#
# Distributed under terms of the GNU GPLv3 license.
# 

import h5py
import numpy
import sys

h5 = h5py.File(sys.argv[1], 'r')

box = numpy.array(h5['/particles/atoms/box/edges/value'][-1])
print box

#cl = h5['/connectivity/fpl/value'][-1]

bonds = [(b1, b2) for cl in h5['/connectivity/fpl/value'] for b1, b2 in cl if b1 != -1 and b2 != -1]

pos = h5['/particles/atoms/position/value'][-1]
b_pos = [(pos[b1], pos[b2]) for b1, b2 in bonds]

def calc_distance(p1, p2, box, half_box):
    d = p2 - p1
    for i in range(3):
        if d[i] > half_box[i]:
            d[i] -= box[i]
        elif d[i] < -half_box[i]:
            d[i] += box[i]
    return d.dot(d)

half_box = 0.5*box
dist = [calc_distance(p1, p2, box, half_box) for p1, p2 in b_pos]

from matplotlib import pyplot as plt

plt.hist(dist, bins=50)
plt.show()
