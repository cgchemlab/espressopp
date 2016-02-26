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

import numpy as np

import logging
import random
import time

import tools
import conf


def _args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--debug', default=False, action='store_true')
    parser.add_argument('--steps', default=100, type=int)
    parser.add_argument('--loops', default=10, type=int)
    parser.add_argument('--rate', default=1.0, type=float)
    parser.add_argument('--interval', default=100, type=int)
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
    args = _args().parse_args()

    system, integrator = tools.initialize_system(args, conf)

    # Prepares system.
    bonds_a_c_tmp, particle_ids = tools.prepare_system(
        conf, system, active_sites=conf.active_sites)

    print('Setup interactions...')
    dynamic_ex_list = espressopp.DynamicExcludeList(integrator, bonds_a_c_tmp)
    verletList = espressopp.VerletList(
        system,
        cutoff=conf.rc)
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
    system.addInteraction(interHarmonic, 'harmonic')

    # Bending potential.
    tpl_a_a = espressopp.FixedTripleList(system.storage)
    potAngHarmonic = espressopp.interaction.AngularHarmonic(
        K=conf.kAng,
        theta0=np.deg2rad(conf.angle),
    )
    interAngHarmonic = espressopp.interaction.FixedTripleListAngularHarmonic(
        system,
        tpl_a_a,
        potAngHarmonic)
    system.addInteraction(interAngHarmonic, 'angle')

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
    system.addInteraction(interLJ, 'lj')

    print('Set max VL cutoff to: {}'.format(max(vl_cutoffs)))
    verletList.setVerletCutoff(max(vl_cutoffs))
    print(verletList.getVerletCutoff())

    total_velocity = espressopp.analysis.TotalVelocity(system)
    total_velocity.reset()

    if args.eq_conf:
        tools.load_conf(system, args.eq_conf)
    else:
        integrator.step = 0
        print('Equilibrating...')
        for i in range(args.eq_loops):
            integrator.run(100)
            tools.info(conf, system, integrator)
        # Save position and velocity of equlibrated configuration
        eq_file = '{}_{}_{}_{}_eq.pck'.format(
            args.prefix, args.rate, args.alpha, args.seed)
        print('Saving equilibrated configuration to {}'.format(eq_file))
        tools.save_conf(system, particle_ids, eq_file)

    if conf.active_sites > 0:
        print('Activating {} molecules'.format(conf.active_sites))
        sample_ids = random.sample(particle_ids, conf.active_sites)
        for pid in sample_ids:
            system.storage.modifyParticle(pid, 'state', 3)

    # Define chemical reaction. fpl_a_a stores new bonds.
    print('Setup chemical reactions...')

    logging.getLogger('ChemicalReaction').setLevel(logging.INFO)

    ar = espressopp.integrator.ChemicalReaction(
        system,
        verletList,
        system.storage,
        args.interval)
    # Reaction: A + A -> B:B + C
    r_type_1 = espressopp.integrator.Reaction(
        type_1=conf.type_a.type_id,
        type_2=conf.type_a.type_id,
        delta_1=1,
        delta_2=-1,
        min_state_1=3,
        max_state_1=4,
        min_state_2=2,
        max_state_2=3,
        rate=args.rate,
        fpl=fpl_a_a,
        intramolecular=True,
        cutoff=1.1*conf.type_a.sigma)
    # conf.rc_lj*tools.lb_sigma(conf.type_a.sigma, conf.type_b.sigma))
    print('Adding reaction 1, rate={}, cutoff={}, P={}'.format(args.rate, r_type_1.cutoff,
                                                               args.rate*args.interval*conf.dt))
    r_type_2 = espressopp.integrator.Reaction(
        type_1=conf.type_a.type_id,
        type_2=conf.type_a.type_id,
        delta_1=-1,
        delta_2=-1,
        min_state_1=2,
        max_state_1=3,
        min_state_2=1,
        max_state_2=2,
        rate=args.rate,
        fpl=fpl_a_a,
        intramolecular=True,
        cutoff=1.1*conf.type_a.sigma)
    # Release one particle
    ar.add_reaction(r_type_1)
    ar.add_reaction(r_type_2)

    r_type_3 = espressopp.integrator.DissociationReaction(
        type_1=conf.type_a.type_id,
        type_2=conf.type_a.type_id,
        delta_1=1,
        delta_2=1,
        min_state_1=0,
        max_state_1=1,
        min_state_2=1,
        max_state_2=2,
        rate=args.rate,
        diss_rate=0.5*args.rate,
        fpl=fpl_a_a,
        cutoff=5.0*conf.type_a.sigma)
    ar.add_reaction(r_type_3)

    r_type_4 = espressopp.integrator.DissociationReaction(
        type_1=conf.type_a.type_id,
        type_2=conf.type_a.type_id,
        delta_1=-1,
        delta_2=1,
        min_state_1=4,
        max_state_1=5,
        min_state_2=1,
        max_state_2=2,
        rate=args.rate,
        diss_rate=0.5*args.rate,
        fpl=fpl_a_a,
        cutoff=5.0)
    ar.add_reaction(r_type_4)

    integrator.addExtension(ar)

    # logging.getLogger('ChemicalReaction').setLevel(logging.DEBUG)
    # logging.getLogger('Reaction').setLevel(logging.DEBUG)
    # logging.getLogger('DissociationReaction').setLevel(logging.DEBUG)
    # logging.getLogger('TopologyManager').setLevel(logging.DEBUG)

    # topology_manager = espressopp.integrator.TopologyManager(system)
    # topology_manager.rebuild()
    # topology_manager.observe_tuple(fpl_a_a)
    # topology_manager.initialize_topology()
    # topology_manager.register_triplet(tpl_a_a, conf.type_a.type_id)

    #integrator.addExtension(topology_manager)

    output_file = '{}_{}_{}_{}.h5'.format(args.prefix, args.rate, args.alpha, args.seed)
    print('Trajectory: {}'.format(output_file))
    traj_file = espressopp.io.DumpH5MD(
        system, output_file,
        group_name='atoms',
        static_box=False,
        author='Jakub Krajniak',
        email='jkrajniak@gmail.com',
        store_species=True,
        store_velocity=True,
        store_state=True,
        store_lambda=True)

    energy_file = '{}_{}_{}_{}_energy.csv'.format(
        args.prefix, args.rate, args.alpha, args.seed)
    print('Energy saved to: {}energy.csv'.format(args.prefix))
    system_analysis = espressopp.analysis.SystemMonitor(
        system, integrator, espressopp.analysis.SystemMonitorOutputCSV(
            energy_file))
    temp = espressopp.analysis.Temperature(system)
    for t in conf.type_ids:
        temp.add_type(t)
    system_analysis.add_observable('T', temp)
    system_analysis.add_observable(
        'Ekin', espressopp.analysis.KineticEnergy(system, temp))
    for label, interaction in system.getAllInteractions().items():
        print('System analysis: adding {}'.format(label))
        system_analysis.add_observable(
            label, espressopp.analysis.PotentialEnergy(system, interaction))

    system_analysis.add_observable(
        'fpl', espressopp.analysis.NFixedPairListEntries(system, fpl_a_a))

    ext_analysis = espressopp.integrator.ExtAnalyze(system_analysis, 100)
    integrator.addExtension(ext_analysis)
    print('Configured system analysis')

    # Save parameters of the simulation
    traj_file.set_parameters({
        'temperature': conf.T,
        'thermostat-gamma': conf.gamma,
        'th': args.interval,
        'rate': args.rate,
        'force-cap': conf.force_cap,
        'alpha': args.alpha,
        'ar-cutoff': r_type_1.cutoff,
        'rng-seed': args.seed,
        'steps': args.steps,
        'loops': args.loops
        })

# Observe tuple
    dump_topol = espressopp.io.DumpTopology(system, integrator, traj_file)
    dump_topol.observe_tuple(fpl_a_a, 'fpl')
    dump_topol.dump()
    dump_topol.update()
    ext_dump = espressopp.integrator.ExtAnalyze(dump_topol, 10)
    integrator.addExtension(ext_dump)


# Reset the integrator step.
    integrator.step = 0
# Run system with non-capped potentials, no thermostat, fixed LJ epsilon
    print('Running serious simulation %s %s in T=%s' % (
        args.loops*args.steps, 'steps', conf.T))
    time0 = time.time()
    traj_file.dump(0, 0)
    system_analysis.info()
    system_analysis.dump()
    r_type_1.active = True
    r_type_2.active = True
    r_type_3.active = True # False #True
    r_type_4.active = True # False #True
    for k in range(args.loops):
        integrator.run(args.steps)
        system_analysis.info()
        traj_file.dump(k*args.steps, k*args.steps*conf.dt)
        traj_file.flush()
        dump_topol.update()
        if k == args.loops / 2:
            print('Activated r3 and r4')
            r_type_3.active = False
            r_type_4.active = False
    traj_file.close()
    dump_topol.update()
    espressopp.tools.analyse.final_info(system, integrator, verletList, time0, time.time())
    a_a_bonds = fpl_a_a.getBonds()
    a_a_a_angles = tpl_a_a.getTriples()
    import cPickle
    cPickle.dump([a_a_bonds, a_a_a_angles], open('topol.pck', 'wb'))
    print('Created {} bonds'.format(fpl_a_a.totalSize()))
    print('Created {} angles'.format(tpl_a_a.totalSize()))


if __name__ == '__main__':
    main()
