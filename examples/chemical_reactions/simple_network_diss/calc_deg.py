import h5py
import networkx as nx
import sys

h5 = h5py.File(sys.argv[1], 'r')

g = nx.Graph()
bonds = h5['/connectivity/fpl/value'][-1]

for b1, b2 in bonds:
    if -1 not in [b1, b2]:
        g.add_edge(b1, b2)

import IPython
IPython.embed()
