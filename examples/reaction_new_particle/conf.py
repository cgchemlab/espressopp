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
type_b = AtomType(2, 0.8*mass, pow(7.0/8.0, 1.0/3.0)*sigma, epsilon)
type_c = AtomType(3, 0.2*mass, pow(1.0/8.0, 1.0/3.0)*sigma, epsilon)
type_c_tmp = AtomType(4, 0.2*mass, pow(1.0/8.0, 1.0/3.0)*sigma, epsilon)

skin = 0.16*sigma
rc = 2.5*sigma
rc_lj = pow(2.0, 1.0/6.0)
dt = 0.0025
T = 0.5
gamma = 5.0

# Bond A-B
kF = 30.0
R0 = 0.5*sigma

# Size of system.
rho = 0.8
N_a = 100
L = 10.0 # pow(N_a*type_a.mass/rho, 1.0/3.0)
box = (L, L, L)

# Co-partner C
R_ac = 0.5*type_a.sigma

types = [type_a, type_b]
type_ids = [x.type_id for x in types]

potential_matrix = [
    (types[i1].type_id, types[i2].type_id,
     tools.lb_sigma(types[i1].sigma, types[i2].sigma),
     tools.lb_epsilon(types[i1].epsilon, types[i2].epsilon))
    for i1 in range(len(types))
    for i2 in range(i1, len(types))
    ]

# Chemical reaction setup
ar_interval = 100
ar_rate = 0.8
ar_cutoff = rc_lj

print('kA = {}'.format(ar_rate*dt*ar_interval))
