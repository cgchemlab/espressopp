#!/usr/bin/env python
"""
Copyright (C) 2015 Jakub Krajniak <jkrajniak@gmail.com>

This file is distributed under free software licence:
you can redistribute it and/or modify it under the terms of the
GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""

import espressopp  # pylint:disable=F0401
import argparse

try:
    import MPI  # pylint: disable=F0401
except ImportError:
    from mpi4py import MPI

import numpy as np

import logging
import time

import tools
import conf


def _args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--debug', default=False, action='store_true')
    parser.add_argument('--steps', default=100, type=int)
    parser.add_argument('--loops', default=10, type=int)
    parser.add_argument('--rate', default=1.0, type=float)
    parser.add_argument('--interval', default=10, type=int)
    parser.add_argument('--alpha', default=0.00001, type=float)
    parser.add_argument('--warmup_loops', default=100, type=int)
    parser.add_argument('--eq_loops', default=1000, type=int)
    parser.add_argument('--eq_conf', default=None)
    parser.add_argument('--prefix', default='dump')
    parser.add_argument('--seed', default=12345, type=int)
    return parser


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
    args = _args().parse_args()

    # enable_debug(args, 'ChemicalReaction')
    # enable_debug(args, 'Reaction')

    print('Box={}'.format(conf.box))

    # Initialize the espressopp system
    system = espressopp.System()
    system.kb = 1.0
    system.rng = espressopp.esutil.RNG(args.seed)
    system.bc = espressopp.bc.OrthorhombicBC(system.rng, conf.box)
    system.skin = conf.skin

    def info():
        espressopp.tools.analyse.info(system, integrator, per_atom=True,
                                      valid_types=conf.type_ids+[conf.type_c.type_id])

    nodeGrid = espressopp.tools.decomp.nodeGrid(MPI.COMM_WORLD.size)
    cellGrid = espressopp.tools.decomp.cellGrid(conf.box, nodeGrid, conf.rc, conf.skin)
    print(nodeGrid, cellGrid)
    system.storage = espressopp.storage.DomainDecomposition(system, nodeGrid, cellGrid)

    integrator = espressopp.integrator.VelocityVerlet(system)
    integrator.dt = conf.dt

    print('Setup integrator, T={}, gamma={}'.format(conf.T, conf.gamma))
    thermostat = espressopp.integrator.LangevinThermostat(system)
    thermostat.temperature = conf.T
    thermostat.gamma = conf.gamma
    thermostat.add_valid_types(conf.type_ids)
    thermostat.add_valid_type_id(conf.type_c.type_id)
    integrator.addExtension(thermostat)

    # Prepares system.
    bonds_a_c_tmp, particle_ids = tools.prepare_system(
        conf, system, active_sites=conf.active_sites)

    # Co-partner C is keep on certain distance. Release cause that molecule is activated.
    fix_list = [(x[0], x[1], conf.R_ac) for x in bonds_a_c_tmp]
    fix_distance = espressopp.integrator.FixDistances(system, fix_list)
    pp = espressopp.integrator.PostProcessChangeProperty()
    pp.add_change_property(
        conf.type_c_tmp.type_id,
        espressopp.ParticleProperties(conf.type_c.type_id, conf.type_c.mass, 0.0))
    fix_distance.add_postprocess(pp)

    integrator.addExtension(fix_distance)

    print('Setup interactions...')
    dynamic_ex_list = espressopp.DynamicExcludeList(integrator, bonds_a_c_tmp)
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

    # Bending potential
    tpl_a_a = espressopp.FixedTripleList(system.storage)
    potAngHarmonic = espressopp.interaction.AngularHarmonic(
        K=conf.kAng,
        theta0=np.deg2rad(conf.angle),
    )
    interAngHarmonic = espressopp.interaction.FixedTripleListAngularHarmonic(
        system,
        tpl_a_a,
        potAngHarmonic)
    system.addInteraction(interAngHarmonic)

    dynamic_ex_list.observe(fpl_a_a)

    if not args.eq_conf:
        tools.warmup(system, integrator, verletList, args, conf)

    # Lennard-Jones potential
    vl_cutoffs = []
    interLJ = espressopp.interaction.VerletListLennardJones(verletList)
    for type_1, type_2, sigma_12, epsilon_12 in conf.potential_matrix:
        print('LJ {}-{}'.format(type_1, type_2))
        interLJ.setPotential(
            type1=type_1,
            type2=type_2,
            potential=espressopp.interaction.LennardJones(
                sigma=sigma_12,
                epsilon=epsilon_12,
                cutoff=sigma_12*conf.rc_lj
            ))
        vl_cutoffs.append(sigma_12*conf.rc_lj)
    system.addInteraction(interLJ)

    # LJ force capped for dummy water molecules
    interLJF = espressopp.interaction.VerletListLennardJonesForceCapped(verletList)
    pot_dummy = espressopp.interaction.LennardJonesForceCapped(
        sigma=conf.type_c_tmp.sigma, epsilon=conf.type_c_tmp.epsilon,
        cutoff=conf.type_c_tmp.sigma*conf.rc_lj)
        cutoff=conf.type_c_tmp.sigma*conf.rc_lj)
    pot_dummy.max_force = conf.force_cap
    interLJF.setPotential(type1=conf.type_c_tmp.type_id, type2=conf.type_c_tmp.type_id,
                          potential=pot_dummy)
        type2=conf.type_c_tmp.type_id,
        potential=pot_dummy)
    system.addInteraction(interLJF)

    # DynamicResolution of Lennard-Jones potential
    interDynamicResLJ = espressopp.interaction.VerletListDynamicResolutionLennardJonesForceCapped(
        verletList, False)
    # Adds C-C interaction.
    pot = espressopp.interaction.LennardJonesForceCapped(
        sigma=conf.type_c.sigma,
        epsilon=conf.type_c.epsilon,
        cutoff=conf.type_c.sigma*conf.rc_lj)
    pot.max_force = conf.force_cap
    interDynamicResLJ.setPotential(
        type1=conf.type_c.type_id, type2=conf.type_c.type_id,
        potential=pot)
    vl_cutoffs.append(conf.type_c.sigma*conf.rc_lj)

    # Adds C - type_ interaction
    for type_id, _, sigma, epsilon in conf.types:
        sigma_12 = tools.lb_sigma(conf.type_c.sigma, sigma)
        epsilon_12 = tools.lb_epsilon(conf.type_c.epsilon, epsilon)
        pot = espressopp.interaction.LennardJonesForceCapped(
            sigma=sigma_12, epsilon=epsilon_12, cutoff=sigma_12*conf.rc_lj)
        pot.max_force = conf.force_cap
        interDynamicResLJ.setPotential(
            type1=conf.type_c.type_id,
            type2=type_id,
            potential=pot)
        vl_cutoffs.append(sigma_12*conf.rc_lj)
    system.addInteraction(interDynamicResLJ)

    print('Set max VL cutoff to: {}'.format(max(vl_cutoffs)))
    verletList.setVerletCutoff(max(vl_cutoffs))
    print(verletList.getVerletCutoff())

    total_velocity = espressopp.analysis.TotalVelocity(system)
    total_velocity.reset()

    if args.eq_conf:
        tools.load_conf(system, args.eq_conf)
    else:
        integrator.step = 0
        print('Equilibratin...')
        for i in range(args.eq_loops):
            integrator.run(500)
            info()
        # Save position and velocity of equlibrated configuration
        print('Saving to equlibrated_conf.pck')
        tools.save_conf(system, particle_ids, 'equilibrated_conf.pck')

    # Define chemical reaction. fpl_a_a stores new bonds.
    print('Setup chemical reactions...')

    ar = espressopp.integrator.ChemicalReaction(
        system,
        verletList,
        fpl_a_a,
        system.storage,
        args.interval)
    # Reaction: A + A -> B:B + C
    r_type_1 = espressopp.integrator.Reaction(
        type_1=conf.type_a.type_id,
        type_2=conf.type_a.type_id,
        delta_1=-1,
        delta_2=-1,
        min_state_1=1,
        max_state_1=4,
        min_state_2=1,
        max_state_2=4,
        rate=args.rate,
        cutoff=1.1*conf.type_a.sigma)
    # conf.rc_lj*tools.lb_sigma(conf.type_a.sigma, conf.type_b.sigma))
    print('Adding reaction 1, rate={}, cutoff={}, P={}'.format(args.rate, r_type_1.cutoff,
                                                               args.rate*args.interval*conf.dt))
    # Release one particle
    r_1_release = espressopp.integrator.PostProcessReleaseParticles(fix_distance, 1)
    r_type_1.add_postprocess(r_1_release, 1)

    ar.add_reaction(r_type_1)
    integrator.addExtension(ar)

    # Dynamic resolution
    basic_dynamic_res = espressopp.integrator.BasicDynamicResolution(
        system, {conf.type_c.type_id: args.alpha})
    integrator.addExtension(basic_dynamic_res)
    d_1_pp = espressopp.integrator.PostProcessChangeProperty()
    d_1_pp.add_change_property(
        conf.type_c.type_id,
        espressopp.ParticleProperties(
            conf.type_c_new.type_id, conf.type_c_new.mass, 0.0))
    basic_dynamic_res.add_postprocess(d_1_pp, 1)

    topology_manager = espressopp.integrator.TopologyManager(system)
    topology_manager.rebuild()
    topology_manager.observe_tuple(fpl_a_a)
    topology_manager.initialize_topology()
    topology_manager.register_triplet(tpl_a_a, conf.type_a.type_id)

    integrator.addExtension(topology_manager)

    output_file = '{}_{}_{}_{}.h5'.format(args.prefix, args.rate, args.alpha, args.seed)
    traj_file = chain_h5md.DumpH5MD(
        output_file,
        system,
        integrator,
        'Jakub Krajniak',
        edges=[edge_i for edge_i in system.bc.boxL], n_states=4,
        valid_types=conf.type_ids + [conf.type_c.type_id])
    topo_file = chain_h5md.DumpTopo(
        traj_file, 'atoms', 'chains', system, integrator, fpl_a_a, time=True,
        chunks=(8, 128, 2))

    # Save parameters of the simulation
    traj_file.parameters.attrs['dt'] = conf.dt
    traj_file.parameters.attrs['skin'] = conf.skin
    traj_file.parameters.attrs['temperature'] = conf.T
    traj_file.parameters.attrs['thermostat-gamma'] = conf.gamma
    traj_file.parameters.attrs['th'] = args.interval
    traj_file.parameters.attrs['rate'] = args.rate
    traj_file.parameters.attrs['force-cap'] = conf.force_cap
    traj_file.parameters.attrs['alpha'] = args.alpha
    traj_file.parameters.attrs['ar-cutoff'] = r_type_1.cutoff
    traj_file.parameters.attrs['rng-seed'] = args.seed
    traj_file.parameters.attrs['steps'] = args.steps
    traj_file.parameters.attrs['loops'] = args.loops
    traj_file.parameters.attrs['active-sites'] = conf.active_sites

# Reset the integrator step.
    integrator.step = 0
    info()
# Run system with non-capped potentials, no thermostat, fixed LJ epsilon
    print('Running serious simulation %s %s in T=%s' % (
        args.loops*args.steps, 'steps', conf.T))
    time0 = time.time()
    T_comp = espressopp.analysis.Temperature(system)
    T_comp.add_type(conf.type_c.type_id)
    # logging.getLogger().setLevel(logging.DEBUG)
    # logging.getLogger("ChemicalReaction").setLevel(logging.DEBUG)
    # logging.getLogger('TopologyManager').setLevel(logging.DEBUG)
    traj_file.analyse()
    traj_file.dump()
    topo_file.dump()

    for k in range(args.loops):
        integrator.run(args.steps)
        info()
        traj_file.analyse()
        traj_file.dump()
        topo_file.dump()
        if k % 10 == 0:
            traj_file.flush()

    espressopp.tools.analyse.final_info(system, integrator, verletList, time0, time.time())
    a_a_bonds = fpl_a_a.getBonds()
    print('Created %d bonds' % len(a_a_bonds[0]))


if __name__ == '__main__':
    main()
