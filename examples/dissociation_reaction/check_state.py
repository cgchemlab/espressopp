#! /usr/bin/env python
#
# Copyright (c) 2016 Jakub Krajniak <jkrajniak@gmail.com>
#
# Distributed under terms of the GNU GPLv3 license.
# 

import argparse
import collections
import h5py
import numpy as np

from matplotlib import pyplot as plt

parser = argparse.ArgumentParser()
parser.add_argument('h5')

args = parser.parse_args()

h5 = h5py.File(args.h5, 'r')

st = h5['/particles/atoms/state/value']

st_time = {x: [] for x in range(np.min(st), np.max(st)+1)}
for s in st:
    d = dict(zip(*map(tuple, np.unique(s, return_counts=True))))
    for k in st_time:
        if k in d:
            st_time[k].append(d[k])
        else:
            st_time[k].append(0)

for k, l in st_time.items():
    if k < 0:
        continue
    plt.plot(l, label=k)
plt.legend()
plt.show()
