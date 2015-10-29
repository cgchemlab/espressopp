"""Test case for releasing new particle."""

from __future__ import print_function

import espressopp  # pylint:disable=F0401
import argparse

try:
    import MPI  # pylint: disable=F0401
except ImportError:
    from mpi4py import MPI
# import numpy

import logging
import math
import time

import serial_h5md

import tools


def _args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--debug', default=False, action='store_true')
    parser.add_argument('--dt', default=0.001, type=float)
    parser.add_argument('--steps', default=100, type=int)
    parser.add_argument('--loops', default=10, type=int)
    parser.add_argument('--warmup_loops', default=100, type=int)
    parser.add_argument('--noinfo', action='store_true', default=False)

    return parser.parse_args()


def _lb_rule_sigma(sigma_a, sigma_b):
    return 0.5*(sigma_a + sigma_b)


def _lb_rule_epsilon(epsilon_a, epsilon_b):
    return math.sqrt(epsilon_a * epsilon_b)


def enable_debug(args, logger=None):
    if args.debug:
        if logger is None:
            logger = ""
        logging.getLogger(logger).setLevel(logging.DEBUG)


def disable_debug():
    logging.getLogger("").setLevel(logging.ERROR)


def main():  # NOQA
    # Modeling reaction
    # A:D + B -> C + D
    args = _args()

    enable_debug(args, 'ChemicalReaction')
    enable_debug(args, 'Reaction')

    epsilon = 1.0
    rc = 2.5  # 2.0**(1./6.)
    rc_lj = 2.0**(1.0/6.0)

    type_a = 1
    type_b = 2
    type_c = 3
    type_d = 4
    type_d_tmp = 5

    sigma_a = sigma_b = 1.0
    sigma_c = 1.0
    sigma_d = 0.5

    mass_a = mass_b = 1.0
    mass_c = 1.5
    mass_d = 0.5

    # Bonds
    # Keep the special particle D by stiff bond inside the A bead
    bond_len_ad = 1.0*sigma_a

    # Newly created bond
    K_ac = 0.1 * epsilon
    bond_len_ac = 0.1

    # Number of particles.
    N_a = 50
    N_b = 50

    TotN = N_a*2 + N_b
    rho = 0.85

    temperature = 1.0

    skin = 0.2
    print("skin = %s" % skin)

    L = pow(TotN/rho, 1.0/3.0)
    box = (L, L, L)
    print('Box={}'.format(box))

    # Initialize the espressopp system
    system = tools.System()
    system.kb = 1.0
    system.rng = espressopp.esutil.RNG()
    system.bc = espressopp.bc.OrthorhombicBC(system.rng, box)
    system.skin = skin

    nodeGrid = espressopp.tools.decomp.nodeGrid(MPI.COMM_WORLD.size)
    cellGrid = espressopp.tools.decomp.cellGrid(box, nodeGrid, rc, skin)
    print(nodeGrid, cellGrid)
    system.storage = espressopp.storage.DomainDecomposition(system, nodeGrid, cellGrid)

    integrator = espressopp.integrator.VelocityVerlet(system)
    integrator.dt = args.dt

    # Build the configuration
    particles_list = []
    ex_list = []
    last_pid = 0

    bonds_a_d = []

    vx, vy, vz = espressopp.tools.init_cfg.velocities.gaussian(
        temperature,
        N_a + N_b,
        [mass_a for _ in range(N_a)] + [mass_b for _ in range(N_b)]
    )

    v_idx = 0
    # Adds A with particle with copartner D.
    for i in range(N_a):
        pos = system.bc.getRandomPos()
        pid_a = last_pid
        vel = espressopp.Real3D(vx[v_idx], vy[v_idx], vz[v_idx])
        particles_list.append([
            pid_a,
            pos,
            vel,
            type_a,
            mass_a,
            1,
            pid_a])
        last_pid += 1
        # Create copartner D
        pid_d = last_pid
        pos_d = pos + espressopp.Real3D(bond_len_ad, 0, 0)
        particles_list.append([pid_d, pos_d, vel, type_d_tmp, mass_d, 1, pid_a])
        last_pid += 1
        bonds_a_d.append((pid_a, pid_d))
        v_idx += 1

    # Adds B
    for i in range(N_b):
        pos = system.bc.getRandomPos()
        vel = espressopp.Real3D(vx[v_idx], vy[v_idx], vz[v_idx])
        pid_b = last_pid
        particles_list.append([pid_b, pos, vel, type_b, mass_b, 0, pid_b])
        last_pid += 1
        v_idx += 1

    system.storage.addParticles(
        particles_list,
        'id',
        'pos',
        'v',
        'type',
        'mass',
        'state',
        'res_id'
        )
    print("Decompose...")
    system.storage.decompose()

    ex_list = bonds_a_d[:]

    logging.getLogger('FixDistances').setLevel(logging.DEBUG)
    fix_list = [(x[0], x[1], bond_len_ad) for x in bonds_a_d]
    fix_distance = espressopp.integrator.FixDistances(system, fix_list, type_a, type_d_tmp)
    fx_1_post_process = espressopp.integrator.PostProcessChangesProperty()
    fx_1_post_process.add_change_property(
        type_d_tmp,
        espressopp.ParticleProperties(type_d, mass_d, 0.0))
    fix_distance.add_postprocess(fx_1_post_process)

    integrator.addExtension(fix_distance)

    print('Setup interactions...')
    verletList = espressopp.VerletList(
        system,
        cutoff=rc,
        exclusionlist=ex_list)

    fpl_a_a = espressopp.FixedPairList(system.storage)
    potHarmonic = espressopp.interaction.Harmonic(
        K=K_ac,
        r0=bond_len_ac,
        cutoff=rc
    )
    interHarmonic = espressopp.interaction.FixedPairListHarmonic(
        system,
        fpl_a_a,
        potHarmonic)
    fpl_a_a.addBonds([])
    system.addInteraction(interHarmonic, 'reaction_a_b')

    # Setup h5md
    h5dump = serial_h5md.DumpH5MD('lj.h5', system, integrator, '', edges=box)
    h5topo = serial_h5md.DumpTopo(
        h5dump,
        system,
        integrator,
        fpl_a_a,
        [x[0] for x in particles_list],
        time=True)

    # Dump existing bonds
    # serial_h5md.save_bonds(h5dump, bond_lists_a + bond_list)
    h5dump.dump()
    h5topo.dump()
    h5dump.f.flush()

    # equilibration
    print('Equilibration...')

    integrator.step = 0
    interEqLJ = espressopp.interaction.VerletListLennardJones(verletList)
    system.addInteraction(interEqLJ, 'lj-warmup')
    eps_delta = 0.0001
    eq_delta = [
        (sigma - eps_delta)/args.warmup_loops
        for _, _, sigma, _ in lj_potential_matrix
        ]
    for s in range(args.warmup_loops):
        espressopp.tools.analyse.info(system, integrator, per_atom=True)
        for i, (type_1, type_2, sigma_12, epsilon_12) in enumerate(lj_potential_matrix):
            interEqLJ.setPotential(
                type1=type_1,
                type2=type_2,
                potential=espressopp.interaction.LennardJones(
                    sigma=eps_delta + s*eq_delta[i],
                    epsilon=epsilon_12,
                    cutoff=(eps_delta+s*eq_delta[i])*rc_lj
                ))
        integrator.run(100)

    system.removeInteractionByLabel('lj-warmup')
    print('Finished equilibration...')

    # Lennard-Jones potential
    interLJ = espressopp.interaction.VerletListLennardJones(verletList)
    for type_1, type_2, sigma_12, epsilon_12 in lj_potential_matrix:
        interLJ.setPotential(
            type1=type_1,
            type2=type_2,
            potential=espressopp.interaction.LennardJones(
                sigma=sigma_12,
                epsilon=epsilon_12,
                cutoff=sigma_12*rc_lj
            ))
    system.addInteraction(interLJ, 'lj')

    interLJNon = espressopp.interaction.VerletListNonReciprocalLennardJones(verletList, type_d_tmp)
    for type_1, type_2, sigma_12, epsilon_12 in lj_potential_non_reciprocal:
        interLJNon.setPotential(
            type1=type_1,
            type2=type_2,
            potential=espressopp.interaction.LennardJones(
                sigma=sigma_12,
                epsilon=epsilon_12,
                cutoff=sigma_12*rc_lj
            ))
    system.addInteraction(interLJNon, 'lj-non')


    print('Equilibrating for 10000')
    for i in range(100):
        integrator.run(1000)
        h5dump.dump()


    # Define chemical reaction. fpl_a_a stores new bonds.
    print('Setup chemical reactions...')

    ar_interval = 100
    ar_rate = 5.0
    ar_cutoff = 1.2*rc_lj
    print('AR interval: {}, AR rate: {}, AR cutoff: {}'.format(
        ar_interval, ar_rate, ar_cutoff))
    ar = espressopp.integrator.ChemicalReaction(
        system,
        verletList,
        fpl_a_a,
        system.storage,
        ar_interval
        )
    r_type_1 = espressopp.integrator.Reaction(
        type_1=type_a,
        type_2=type_b,
        delta_1=0,
        delta_2=1,
        min_state_1=1,
        max_state_1=2,
        min_state_2=0,
        max_state_2=1,
        rate=ar_rate,
        cutoff=ar_cutoff
    )
    r_1_post_process = espressopp.integrator.PostProcessChangesProperty()
    r_1_post_process.add_change_property(
        type_a,
        espressopp.ParticleProperties(type_c, mass_c, 0.0))
    r_type_1.add_postprocess(r_1_post_process)

    r_type_2 = espressopp.integrator.Reaction(
        type_1=type_b,
        type_2=type_a,
        delta_1=-2,
        delta_2=0,
        min_state_1=1,
        max_state_1=2,
        min_state_2=1,
        max_state_2=2,
        rate=ar_rate,
        cutoff=ar_cutoff
    )
    r_2_post_process = espressopp.integrator.PostProcessChangesProperty()
    r_2_post_process.add_change_property(
        type_a,
        espressopp.ParticleProperties(type_c, mass_c, 0.0))
    r_type_2.add_postprocess(r_2_post_process)

    r_type_3 = espressopp.integrator.Reaction(
        type_1=type_c,
        type_2=type_b,
        delta_1=-1,
        delta_2=1,
        min_state_1=1,
        max_state_1=2,
        min_state_2=0,
        max_state_2=1,
        rate=ar_rate,
        cutoff=ar_cutoff
    )

    r_type_4 = espressopp.integrator.Reaction(
        type_1=type_b,
        type_2=type_c,
        delta_1=-2,
        delta_2=-1,
        min_state_1=1,
        max_state_1=2,
        min_state_2=1,
        max_state_2=2,
        rate=ar_rate,
        cutoff=ar_cutoff
    )
    ar.add_reaction(r_type_1)
    ar.add_reaction(r_type_2)
    ar.add_reaction(r_type_3)
    ar.add_reaction(r_type_4)

    integrator.addExtension(ar)
    print('Chemical reaction with the rate={}, dt={}, interval={}'.format(
        ar_rate, args.dt, ar_interval))

    thermostat = espressopp.integrator.StochasticVelocityRescaling(system)
    thermostat.temperature = temperature
    thermostat.coupling = 0.5
    integrator.addExtension(thermostat)

    print('Potential description')
    system.printInteractions()

# Reset the integrator step.
    integrator.step = 0
    espressopp.tools.analyse.info(system, integrator, per_atom=True)
# Run system with non-capped potentials, no thermostat, fixed LJ epsilon
    print('Running serious simulation %s %s in T=%s' % (
        args.loops*args.steps, 'steps', temperature))
    time0 = time.time()
    import sys
    for k in range(args.loops):
        integrator.run(args.steps)
        if not args.noinfo:
            sys.stdout.write('[ {} ] '.format(len(fpl_a_a.getBonds()[0])))
            espressopp.tools.analyse.info(system, integrator, per_atom=True)
        h5dump.dump()
        h5topo.dump()
        h5dump.analyse()
        h5dump.f.flush()

    h5dump.close()
    espressopp.tools.analyse.final_info(system, integrator, verletList, time0, time.time())
    a_a_bonds = fpl_a_a.getBonds()
    print('Created %d bonds' % len(a_a_bonds[0]))
    print('Bonds:')
    print('\n'.join(map(str, a_a_bonds[0])))

    print('Index')
    print('index {}'.format(' '.join(map(str, set([p for x in a_a_bonds[0] for p in x])))))


if __name__ == '__main__':
    main()
