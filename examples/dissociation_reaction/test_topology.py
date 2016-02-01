#! /usr/bin/env python
#
# Copyright (c) 2016 Jakub Krajniak <jkrajniak@gmail.com>
#
# Distributed under terms of the GNU GPLv3 license.
# 

import networkx as nx
import cPickle
import sys

from matplotlib import pyplot as plt

topol_data = cPickle.load(open(sys.argv[1], 'rb'))
print len(topol_data[0]), len(topol_data[1])

bonds, angles, bonds_1_4 = topol_data

print angles
print bonds_1_4

g = nx.Graph()
for bl in bonds:
    for b1, b2 in bl:
        g.add_edge(b1, b2)

nodes = g.nodes()
gen_angles = set()
for idx, i in enumerate(nodes):
    for j in nodes[idx+1:]:
        path_i_j = nx.all_simple_paths(g, i, j, 4)
        for x in path_i_j:
            if len(x) == 3:
                a, ra = tuple(x), tuple(reversed(x))
                gen_angles.add(a)
                gen_angles.add(ra)

u_angles = set()
eu_angles = set()
for al in angles:
    for a in al:
        u_angles.add(a)
        eu_angles.add(a)
        eu_angles.add(tuple(reversed(a)))

print len(gen_angles)
print len(eu_angles)

print gen_angles - eu_angles
