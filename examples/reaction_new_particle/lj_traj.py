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

import tools
import conf


def _args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--debug', default=False, action='store_true')
    parser.add_argument('--steps', default=100, type=int)
    parser.add_argument('--loops', default=10, type=int)
    parser.add_argument('--warmup_loops', default=100, type=int)
    parser.add_argument('--eq_conf', default=None)
    parser.add_argument('--vis', action='store_true', default=False)

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

    # enable_debug(args, 'ChemicalReaction')
    # enable_debug(args, 'Reaction')

    print('Box={}'.format(conf.box))

    # Initialize the espressopp system
    system = espressopp.System()
    system.kb = 1.0
    system.rng = espressopp.esutil.RNG()
    system.bc = espressopp.bc.OrthorhombicBC(system.rng, conf.box)
    system.skin = conf.skin

    def info():
        espressopp.tools.analyse.info(system, integrator, per_atom=True,
                                      valid_types=conf.type_ids)

    nodeGrid = espressopp.tools.decomp.nodeGrid(MPI.COMM_WORLD.size)
    cellGrid = espressopp.tools.decomp.cellGrid(conf.box, nodeGrid, conf.rc, conf.skin)
    print(nodeGrid, cellGrid)
    system.storage = espressopp.storage.DomainDecomposition(system, nodeGrid, cellGrid)

    integrator = espressopp.integrator.VelocityVerlet(system)
    integrator.dt = conf.dt

    bonds_a_c_tmp, particle_ids = tools.prepare_system(conf, system)

    ex_list = bonds_a_c_tmp[:]
    # Co-partner C is keep on certain distance.
    fix_list = [(x[0], x[1], conf.R_ac) for x in bonds_a_c_tmp]
    fix_distance = espressopp.integrator.FixDistances(
        system, fix_list, conf.type_a.type_id, conf.type_c_tmp.type_id)
    fx_1_post_process = espressopp.integrator.PostProcessChangeProperty()
    fx_1_post_process.add_change_property(
        conf.type_c_tmp.type_id,
        espressopp.ParticleProperties(conf.type_c.type_id, conf.type_c.mass, 0.0))
    fix_distance.add_postprocess(fx_1_post_process)

    integrator.addExtension(fix_distance)

    print('Setup interactions...')
    dynamic_ex_list = espressopp.DynamicExcludeList(integrator, ex_list)
    verletList = espressopp.VerletList(
        system,
        cutoff=conf.rc,
        exclusionlist=dynamic_ex_list)

    # Chemical bond A-B and B-B.
    fpl_a_a = espressopp.FixedPairList(system.storage)
    potHarmonic = espressopp.interaction.Harmonic(
        K=conf.kF,
        r0=conf.R0,
        cutoff=conf.rc)
    interHarmonic = espressopp.interaction.FixedPairListHarmonic(
        system,
        fpl_a_a,
        potHarmonic)
    fpl_a_a.addBonds([])
    system.addInteraction(interHarmonic)

    if not args.eq_conf:
        tools.warmup(system, integrator, verletList, args, conf)

    # Lennard-Jones potential
    interLJ = espressopp.interaction.VerletListLennardJones(verletList)
    for type_1, type_2, sigma_12, epsilon_12 in conf.potential_matrix:
        interLJ.setPotential(
            type1=type_1,
            type2=type_2,
            potential=espressopp.interaction.LennardJones(
                sigma=sigma_12,
                epsilon=epsilon_12,
                cutoff=sigma_12*conf.rc_lj
            ))
    system.addInteraction(interLJ)

    # DynamicResolution of Lennard-Jones potential
    interDynamicResLJ = espressopp.interaction.VerletListDynamicResolutionLennardJones(
        verletList, False)
    interDynamicResLJ.setPotential(
        type1=conf.type_c.type_id, type2=conf.type_c.type_id,
        potential=espressopp.interaction.LennardJones(
            sigma=conf.type_c.sigma,
            epsilon=conf.type_c.epsilon,
            cutoff=conf.sigma*conf.rc_lj))
    for type_id, _, sigma, epsilon in conf.types:
        interDynamicResLJ.setPotential(
            type1=conf.type_c.type_id,
            type2=type_id,
            potential=espressopp.interaction.LennardJones(
                sigma=tools.lb_sigma(conf.type_c.sigma, sigma),
                epsilon=tools.lb_epsilon(conf.type_c.epsilon, epsilon),
                cutoff=conf.sigma*conf.rc_lj))
    system.addInteraction(interDynamicResLJ)

    total_velocity = espressopp.analysis.TotalVelocity(system)
    total_velocity.reset()

    if args.eq_conf:
        tools.load_conf(system, args.eq_conf)
    else:
        print('Equilibratin...')
        for i in range(args.warmup_loops):
            integrator.run(100)
            info()
        # Save position and velocity of equlibrated configuration
        print('Saving to equlibrated_conf.pck')
        tools.save_conf(system, particle_ids, 'equilibrated_conf.pck')

    # Define chemical reaction. fpl_a_a stores new bonds.
    print('Setup chemical reactions...')

    print('AR interval: {}, AR rate: {}, AR cutoff: {}'.format(
        conf.ar_interval, conf.ar_rate, conf.ar_cutoff))
    ar = espressopp.integrator.ChemicalReaction(
        system,
        verletList,
        fpl_a_a,
        system.storage,
        conf.ar_interval)
    # Reaction: A + A -> B:B + C
    r_type_1 = espressopp.integrator.Reaction(
        type_1=conf.type_a.type_id,
        type_2=conf.type_a.type_id,
        delta_1=-1,
        delta_2=-1,
        min_state_1=1,
        max_state_1=3,
        min_state_2=1,
        max_state_2=3,
        rate=conf.ar_rate,
        cutoff=conf.ar_cutoff)
    # Change type: A -> B
    r_1_post_process = espressopp.integrator.PostProcessChangeProperty()
    r_1_post_process.add_change_property(
        conf.type_a.type_id,
        espressopp.ParticleProperties(conf.type_b.type_id, conf.type_b.mass, 0.0))
    r_type_1.add_postprocess(r_1_post_process)

    # Reaction: A + B:B -> B:B:B + C
    r_type_2 = espressopp.integrator.Reaction(
        type_1=conf.type_a.type_id,
        type_2=conf.type_b.type_id,
        delta_1=-1,
        delta_2=-1,
        min_state_1=1,
        max_state_1=3,
        min_state_2=1,
        max_state_2=2,
        rate=conf.ar_rate,
        cutoff=conf.ar_cutoff)
    # Change type: A -> B
    r_2_post_process = espressopp.integrator.PostProcessChangeProperty()
    r_2_post_process.add_change_property(
        conf.type_a.type_id,
        espressopp.ParticleProperties(conf.type_b.type_id, conf.type_b.mass, 0.0))
    r_type_2.add_postprocess(r_2_post_process)

    # Reaction: B + B:B -> B:B:B
    r_type_3 = espressopp.integrator.Reaction(
        type_1=conf.type_b.type_id,
        type_2=conf.type_b.type_id,
        delta_1=-1,
        delta_2=-1,
        min_state_1=1,
        max_state_1=2,
        min_state_2=1,
        max_state_2=2,
        rate=conf.ar_rate,
        cutoff=conf.ar_cutoff)

    ar.add_reaction(r_type_1)
    ar.add_reaction(r_type_2)
    ar.add_reaction(r_type_3)

    integrator.addExtension(ar)
    print('Chemical reaction with the rate={}, dt={}, interval={}'.format(
        conf.ar_rate, conf.dt, conf.ar_interval))

    # Dynamic resolution
    basic_dynamic_res = espressopp.integrator.BasicDynamicResolution(
        system, {conf.type_c.type_id: 0.0001})
    integrator.addExtension(basic_dynamic_res)

    thermostat = espressopp.integrator.StochasticVelocityRescaling(system)
    thermostat.temperature = conf.T
    thermostat.coupling = 0.5
    integrator.addExtension(thermostat)

    dump_gro = espressopp.io.DumpGRO(
        system, integrator, filename='trajectory.gro', append=True
        )
    ext_dump_gro = espressopp.integrator.ExtAnalyze(dump_gro, 10)
    integrator.addExtension(ext_dump_gro)

    ps = espressopp.analysis.PyStore(
        system, 'trj.h5', group_name='atoms', store_state=True,
        store_lambda=True,
        static_box=True, author='Jakub Krajniak', email='jkrajniak@gmail.com')

    ps.add_connectivity(fpl_a_a, 'a_a')

    ps.dump(0, 0)

    # ext_total_vel = espressopp.integrator.ExtAnalyze(total_velocity, 10)
    # integrator.addExtension(ext_total_vel)

    if args.vis:
        import networkx as nx
        from matplotlib import pyplot as plt
        g = nx.Graph()
        nx.draw_graphviz(g)
        plt.show(False)
        plt.draw()

# Reset the integrator step.
    integrator.step = 0
    info()
# Run system with non-capped potentials, no thermostat, fixed LJ epsilon
    print('Running serious simulation %s %s in T=%s' % (
        args.loops*args.steps, 'steps', conf.T))
    time0 = time.time()
    import sys
    print [system.storage.getParticle(x).lambda_adr for x in particle_ids].count(0.0)
    for k in range(args.loops):
        integrator.run(args.steps)
        print [system.storage.getParticle(x).lambda_adr for x in particle_ids].count(0.0)
        sys.stdout.write('[ {} ] '.format(len(fpl_a_a.getBonds()[0])))
        sys.stdout.write('~ {} ~ '.format(fix_distance.size))
        sys.stdout.write('+ {} + '.format(dynamic_ex_list.size))
        info()
        ps.dump(k*args.steps, k*args.steps*conf.dt)

        if args.vis:
            aa_bonds = fpl_a_a.getBonds()[0]
            g.add_edges_from(aa_bonds)
            nx.draw_graphviz(g)
            plt.draw()

    ps.dump(k*args.steps, k*args.steps*conf.dt)
    ps.close()
    espressopp.tools.analyse.final_info(system, integrator, verletList, time0, time.time())
    a_a_bonds = fpl_a_a.getBonds()
    print('Created %d bonds' % len(a_a_bonds[0]))


if __name__ == '__main__':
    main()
