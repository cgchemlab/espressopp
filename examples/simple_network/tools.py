#! /usr/bin/env python
#
# Copyright (c) 2015 Jakub Krajniak <jkrajniak@gmail.com>
#
# Distributed under terms of the GNU GPLv3 license.
#

import cPickle
import random
import espressopp


def lb_epsilon(eps1, eps2):
    return pow(eps1*eps2, 1.0/2.0)


def lb_sigma(sig1, sig2):
    return 0.5*(sig1+sig2)


def save_conf(system, particle_ids, output_file):
    output_list = []
    for x in particle_ids:
        p = system.storage.getParticle(x)
        output_list.append({
            'pid': x, 'pos': p.pos, 'v': p.v,
            'type': p.type, 'mass': p.mass,
            'state': p.state,
            'res_id': p.res_id,
            'lambda_adr': p.lambda_adr
        })
    with open(output_file, 'wb') as f:
        cPickle.dump(output_list, f)


def load_conf(system, input_file):
    with open(input_file, 'rb') as f:
        input_list = cPickle.load(f)
        for l in input_list:
            pid = l['pid']
            for prop, val in l.iteritems():
                if prop != 'pid':
                    system.storage.modifyParticle(pid, prop, val)


def prepare_system(conf, system, active_sites=1):
    # Build the configuration.
    particles_list = []
    last_pid = 0

    bonds_a_c_tmp = []

    vx, vy, vz = espressopp.tools.init_cfg.velocities.gaussian(
        conf.T,
        conf.N_a,
        [conf.type_a.mass for _ in range(conf.N_a)]
    )

    v_idx = 0
    pids = []
    # Adds M molecules
    for i in range(conf.N_a):
        pos = system.bc.getRandomPos()
        pid_a = last_pid
        vel = espressopp.Real3D(vx[v_idx], vy[v_idx], vz[v_idx])
        particles_list.append([
            pid_a,
            pos,
            vel,
            conf.type_a.type_id,
            conf.type_a.mass,
            3,
            pid_a,
            1.0])
        pids.append(i)
        last_pid += 1
        # Create three co-partner C
        pid_c = last_pid
        for _ in range(3):
            pos_c = pos + conf.R_ac*system.rng.uniformOnSphere()
            particles_list.append([
                pid_c,
                pos_c,
                vel,
                conf.type_c_tmp.type_id,
                conf.type_c_tmp.mass,
                4,
                pid_c,
                10**-6])
            last_pid += 1
            bonds_a_c_tmp.append((pid_a, pid_c))
        v_idx += 1
    part_prop = ['id', 'pos', 'v', 'type', 'mass', 'state', 'res_id', 'lambda_adr']
    system.storage.addParticles(particles_list, *part_prop)
    print("Decompose...")
    system.storage.decompose()
    return bonds_a_c_tmp, [x[0] for x in particles_list]


def warmup(system, integrator, verletList, args, conf):
    # Equilibration.
    print('Equilibration...')
    integrator.step = 0
    interEqLJ = espressopp.interaction.VerletListLennardJones(verletList)
    system.addInteraction(interEqLJ)
    eps_delta = 0.0001
    eq_delta = [
        (sigma - eps_delta)/args.warmup_loops
        for _, _, sigma, _ in conf.warmup_potential_matrix
        ]
    potential_matrix = {
        (type_1, type_2): espressopp.interaction.LennardJones()
        for type_1, type_2, _, _ in conf.warmup_potential_matrix
    }

    espressopp.tools.analyse.info(system, integrator, per_atom=True, valid_types=conf.type_ids)
    for s in range(args.warmup_loops):
        integrator.run(100)
        espressopp.tools.analyse.info(system, integrator, per_atom=True, valid_types=conf.type_ids)
        for i, (type_1, type_2, sigma_12, epsilon_12) in enumerate(conf.warmup_potential_matrix):
            sigma = eps_delta + s*eq_delta[i]
            cutoff = sigma * conf.rc_lj
            pot = potential_matrix[(type_1, type_2)]
            pot.sigma = sigma
            pot.epsilon = epsilon_12
            pot.cutoff = cutoff
            interEqLJ.setPotential(
                type1=type_1,
                type2=type_2,
                potential=pot)

    print('Finished warming up')
    system.removeInteraction(1)
