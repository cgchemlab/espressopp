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
type_b = AtomType(2, 0.8*mass, pow(124/125.0, 1.0/3.0)*sigma, epsilon)
# New molecule
type_c = AtomType(3, 0.2*mass, 0.2*sigma, epsilon)
type_c_tmp = AtomType(4, type_c.mass, type_c.sigma, type_c.epsilon)
type_c_final = AtomType(5, type_c.mass, type_c.sigma, type_c.epsilon)

skin = 0.3*sigma
rc = 1.5*sigma
rc_lj = pow(2.0, 1.0/6.0)
dt = 0.0025
T = 0.5
gamma = 5.0

# Bond A-B
kF = 30.0
R0 = 0.8*sigma

# Angle
kAng = 25.0
angle = 120.0

# Size of system.
rho = 0.8
N_a = 250
# Number of coopartners.
N_c = 3
active_sites = 1
L = 10.0
#L = pow(N_a*type_a.mass/rho, 1.0/3.0)
box = (L, L, L)

force_cap = 1000.0

# Co-partner, fixed distance A-C_tmp
R_ac = tools.lb_sigma(type_a.sigma, type_c.sigma)*rc_lj+0.01*sigma
print('R_ac={}'.format(R_ac))

types = [type_a, type_b, type_c_final]
type_ids = [x.type_id for x in types]

potential_matrix = [
    (types[i1].type_id, types[i2].type_id,
     tools.lb_sigma(types[i1].sigma, types[i2].sigma),
     tools.lb_epsilon(types[i1].epsilon, types[i2].epsilon))
    for i1 in range(len(types))
    for i2 in range(i1, len(types))
    ]

warmup_potential_matrix = potential_matrix[:]
warmup_potential_matrix.extend([
    (type_c.type_id, types[i].type_id,
     tools.lb_sigma(type_c.sigma, types[i].sigma),
     tools.lb_epsilon(type_c.epsilon, types[i].epsilon))
    for i in range(len(types))])
warmup_potential_matrix.append([type_c.type_id, type_c.type_id, type_c.sigma, type_c.epsilon])
warmup_potential_matrix.append([type_c_tmp.type_id, type_c_tmp.type_id,
                                type_c_tmp.sigma, type_c_tmp.epsilon])
