#! /usr/bin/env python
#
# Copyright (c) 2015 Jakub Krajniak <jkrajniak@gmail.com>
#
# Distributed under terms of the GNU GPLv3 license.
#

from collections import namedtuple

import tools

AtomType = namedtuple('AtomType', ['type_id', 'mass', 'sigma', 'epsilon'])

epsilon = 1.0
sigma = 1.0
mass = 1.0

type_a = AtomType(1, 1.0*mass, 1.0*sigma, epsilon)

energy_groups = [type_a]

skin = 0.3*sigma
rc = 1.5*sigma
rc_lj = pow(2.0, 1.0/6.0)
dt = 0.001
T = 0.5
gamma = 5.0

# Bond A-B
kF = 30.0
R0 = 0.97*sigma

# Angle
kAng = 25.0
angle = 120.0

# Size of system.
rho = 0.8
N_a = 1000
# Number of coopartners.
N_c = 0
active_sites = int(0.2*N_a)
L = pow(N_a*type_a.mass/rho, 1.0/3.0)
box = (L, L, L)

force_cap = 1000.0

types = [type_a]
type_ids = [x.type_id for x in types]

potential_matrix = [
    (types[i1].type_id, types[i2].type_id,
     tools.lb_sigma(types[i1].sigma, types[i2].sigma),
     tools.lb_epsilon(types[i1].epsilon, types[i2].epsilon))
    for i1 in range(len(types))
    for i2 in range(i1, len(types))
    ]

warmup_potential_matrix = potential_matrix[:]
