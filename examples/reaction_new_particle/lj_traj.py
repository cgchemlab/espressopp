"""Test case for releasing new particle."""

import espressopp  # pylint:disable=F0401
import argparse

try:
    import MPI  # pylint: disable=F0401
except ImportError:
    from mpi4py import MPI
# import numpy

import logging
import time

import conf


def _args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--debug', default=False, action='store_true')
    parser.add_argument('--dt', default=0.001, type=float)
    parser.add_argument('--steps', default=100, type=int)
    parser.add_argument('--loops', default=10, type=int)
    parser.add_argument('--warmup_loops', default=100, type=int)
    parser.add_argument('--noinfo', action='store_true', default=False)

    return parser.parse_args()


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

    print('Box={}'.format(conf.box))

    # Initialize the espressopp system
    system = espressopp.System()
    system.kb = 1.0
    system.rng = espressopp.esutil.RNG()
    system.bc = espressopp.bc.OrthorhombicBC(system.rng, conf.box)
    system.skin = conf.skin

    nodeGrid = espressopp.tools.decomp.nodeGrid(MPI.COMM_WORLD.size)
    cellGrid = espressopp.tools.decomp.cellGrid(conf.box, nodeGrid, conf.rc, conf.skin)
    print(nodeGrid, cellGrid)
    system.storage = espressopp.storage.DomainDecomposition(system, nodeGrid, cellGrid)

    integrator = espressopp.integrator.VelocityVerlet(system)
    integrator.dt = args.dt

    # Build the configuration.
    particles_list = []
    last_pid = 0

    bonds_a_d_tmp = []

    vx, vy, vz = espressopp.tools.init_cfg.velocities.gaussian(
        conf.temperature,
        conf.N_a + conf.N_b,
        [conf.type_a.mass for _ in range(conf.N_a)] + [conf.type_b.mass for _ in range(conf.N_b)]
    )

    v_idx = 0
    # Adds A with particle with copartner D.
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
            1,
            pid_a])
        last_pid += 1
        # Create co-partner D
        pid_d = last_pid
        pos_d = pos + espressopp.Real3D(conf.bond_a_d_tmp, 0, 0)
        particles_list.append([
            pid_d, pos_d, vel, conf.type_d_tmp.type_id,
            conf.type_d_tmp.mass, 1, pid_d])
        last_pid += 1
        bonds_a_d_tmp.append((pid_a, pid_d))
        v_idx += 1

    # Adds B
    for i in range(conf.N_b):
        pos = system.bc.getRandomPos()
        vel = espressopp.Real3D(vx[v_idx], vy[v_idx], vz[v_idx])
        pid_b = last_pid
        particles_list.append([pid_b, pos, vel, conf.type_b.type_id, conf.type_b.mass, 0, pid_b])
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
        'res_id')
    print("Decompose...")
    system.storage.decompose()

    ex_list = bonds_a_d_tmp[:]

    logging.getLogger('FixDistances').setLevel(logging.DEBUG)
    fix_list = [(x[0], x[1], conf.bond_a_d_tmp) for x in bonds_a_d_tmp]
    fix_distance = espressopp.integrator.FixDistances(
        system, fix_list, conf.type_a.type_id, conf.type_d_tmp.type_id)
    fx_1_post_process = espressopp.integrator.PostProcessChangeProperty()
    fx_1_post_process.add_change_property(
        conf.type_d_tmp.type_id,
        espressopp.ParticleProperties(conf.type_d.type_id, conf.type_d.mass, 0.0))
    fix_distance.add_postprocess(fx_1_post_process)

    integrator.addExtension(fix_distance)

    print('Setup interactions...')
    dynamic_ex_list = espressopp.DynamicExcludeList(integrator, ex_list)
    verletList = espressopp.VerletList(
        system,
        cutoff=conf.rc,
        exclusionlist=dynamic_ex_list)

    # Chemical bond C = A-A.
    fpl_a_a = espressopp.FixedPairList(system.storage)
    potHarmonic = espressopp.interaction.Harmonic(
        K=conf.K_ac,
        r0=conf.bond_a_c,
        cutoff=conf.rc
    )
    interHarmonic = espressopp.interaction.FixedPairListHarmonic(
        system,
        fpl_a_a,
        potHarmonic)
    fpl_a_a.addBonds([])
    system.addInteraction(interHarmonic)

    # Equilibration.
    print('Equilibration...')
    integrator.step = 0
    interEqLJ = espressopp.interaction.VerletListLennardJones(verletList)
    system.addInteraction(interEqLJ)
    eps_delta = 0.0001
    eq_delta = [
        (sigma - eps_delta)/args.warmup_loops
        for _, _, sigma, _ in conf.lj_potential_matrix
        ]
    for s in range(args.warmup_loops):
        espressopp.tools.analyse.info(system, integrator, per_atom=True)
        for i, (type_1, type_2, sigma_12, epsilon_12) in enumerate(conf.lj_potential_matrix):
            interEqLJ.setPotential(
                type1=type_1,
                type2=type_2,
                potential=espressopp.interaction.LennardJones(
                    sigma=eps_delta + s*eq_delta[i],
                    epsilon=epsilon_12,
                    cutoff=(eps_delta+s*eq_delta[i])*conf.rc_lj
                ))
        integrator.run(100)

    system.removeInteraction(1)
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
