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
    parser.add_argument('--alpha', default=0.001, type=float)
    parser.add_argument('--warmup_loops', default=100, type=int)
    parser.add_argument('--eq_loops', default=500, type=int)
    parser.add_argument('--eq_conf', default=None)
    parser.add_argument('--prefix', default='dump')
    parser.add_argument('--seed', default=12345, type=int)
    return parser


def main():  # NOQA
    args = _args().parse_args()

    system, integrator = tools.initialize_system(args, conf)

    # Prepares system.
    bonds_a_c_tmp = []
    particles_list = [
        (1, espressopp.Real3D(2.0, 2.0, 2.0), espressopp.Real3D(0.0, 1.0, 0.0),
         conf.type_a.type_id, conf.type_a.mass, 1, 1, 1.0),
        (2, espressopp.Real3D(2.0, 2.0+conf.R0, 2.0),
         espressopp.Real3D(0.0, -1.0, 0.0),
         conf.type_a.type_id, conf.type_a.mass, 2, 1, 1.0),
        (3, espressopp.Real3D(0.1, 2.0+0.1*conf.R0, 2.0),
         espressopp.Real3D(1.5, 0.0, 0.0),
         conf.type_c_final.type_id,
         conf.type_c_final.mass, 4, 2, 1.0),
        (4, espressopp.Real3D(1.41, 2.2, 2.0), espressopp.Real3D(0.0, 1.0, 0.0),
         conf.type_a.type_id, conf.type_a.mass, 2, 1, 1.0)
    ]
    part_prop = ['id', 'pos', 'v', 'type', 'mass', 'state', 'res_id', 'lambda_adr']
    system.storage.addParticles(particles_list, *part_prop)

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
    dynamic_ex_list = espressopp.DynamicExcludeList(integrator, [(1, 2), (1, 4), (2, 4)])
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
    system.addInteraction(interHarmonic, 'harmonic')

    # Exclude those particles from VL after bond is created.
    dynamic_ex_list.observe_tuple(fpl_a_a)
    fpl_a_a.addBonds([(1, 2), (1, 4)])

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
    dynamic_ex_list.observe_triple(tpl_a_a)
    tpl_a_a.addTriples([(1,2,4)])

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

    # LJ force capped for dummy water molecules to prevent those molecules
    # from overlapping.
    interLJF = espressopp.interaction.VerletListSoftCosine(verletList)
    pot_dummy = espressopp.interaction.SoftCosine(cutoff=conf.type_c_tmp.sigma*conf.rc_lj)
    interLJF.setPotential(
        type1=conf.type_c_tmp.type_id,
        type2=conf.type_c_tmp.type_id,
        potential=pot_dummy)
    system.addInteraction(interLJF, 'ljcos-dummy')

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
    system.addInteraction(interDynamicResLJ, 'lj-lmb')

    print('Set max VL cutoff to: {}'.format(max(vl_cutoffs)))
    verletList.setVerletCutoff(max(vl_cutoffs))
    print(verletList.getVerletCutoff())

    # Define chemical reaction. fpl_a_a stores new bonds.
    print('Setup chemical reactions...')

    ar = espressopp.integrator.ChemicalReaction(
        system,
        verletList,
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
        fpl=fpl_a_a,
        cutoff=1.1*conf.type_a.sigma)
    r_type_1.intramolecular = True

    r_type_2 = espressopp.integrator.Reaction(
        type_1=conf.type_c_final.type_id,
        type_2=conf.type_a.type_id,
        delta_1=0,
        delta_2=1,
        min_state_1=4,
        max_state_1=5,
        min_state_2=0,
        max_state_2=3,
        rate=args.rate,
        fpl=fpl_a_a,
        cutoff=0.8)
    r_type_2.intramolecular = True
    r_2_join = espressopp.integrator.PostProcessJoinParticles(fix_distance, conf.R_ac)
    r_type_2.add_postprocess(r_2_join, 'type_1')
    r_2_remove_bond = espressopp.integrator.PostProcessRemoveBond(fpl_a_a, 1)
    r_type_2.add_postprocess(r_2_remove_bond, 'type_1')
    r_2_property = espressopp.integrator.PostProcessChangeProperty()
    r_2_property.add_change_property(
        conf.type_c_final.type_id,
        espressopp.ParticleProperties(
            conf.type_c_tmp.type_id, conf.type_c_final.mass, 0.0, 10**-7))
    r_type_2.add_postprocess(r_2_property, 'type_1')

    # conf.rc_lj*tools.lb_sigma(conf.type_a.sigma, conf.type_b.sigma))
    print('Adding reaction 1, rate={}, cutoff={}, P={}'.format(args.rate, r_type_1.cutoff,
                                                               args.rate*args.interval*conf.dt))
    # Release one particle
    r_1_release = espressopp.integrator.PostProcessReleaseParticles(fix_distance, 1)
    r_type_1.add_postprocess(r_1_release, 1)

    ar.add_reaction(r_type_1)
    ar.add_reaction(r_type_2)
    integrator.addExtension(ar)

    # Dynamic resolution. When particle reaches resolution of 1.0 then
    # it changes type to final type.
    basic_dynamic_res = espressopp.integrator.BasicDynamicResolution(
        system, {conf.type_c.type_id: args.alpha})
    d_1_pp = espressopp.integrator.PostProcessChangeProperty()
    d_1_pp.add_change_property(
        conf.type_c.type_id,
        espressopp.ParticleProperties(
            conf.type_c_final.type_id, conf.type_c_final.mass, 0.0))
    basic_dynamic_res.add_postprocess(d_1_pp, 1)
    integrator.addExtension(basic_dynamic_res)

    topology_manager = espressopp.integrator.TopologyManager(system)
    topology_manager.rebuild()
    topology_manager.observe_tuple(fpl_a_a)
    topology_manager.initialize_topology()
    topology_manager.register_triplet(tpl_a_a, conf.type_a.type_id)
    integrator.addExtension(topology_manager)

    output_file = '{}_{}_{}_{}.h5'.format(args.prefix, args.rate, args.alpha, args.seed)
    print('Trajectory: {}'.format(output_file))
    traj_file = espressopp.io.DumpH5MD(
        system, output_file,
        group_name='atoms',
        static_box=False,
        author='Jakub Krajniak',
        email='jkrajniak@gmail.com',
        store_species=True,
        store_state=True,
        store_res_id=True,
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
        'fpl_a_a', espressopp.analysis.NFixedPairListEntries(system, fpl_a_a))

    ext_analysis = espressopp.integrator.ExtAnalyze(system_analysis, 1)
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

    import logging
    logging.getLogger('Reaction').setLevel(logging.DEBUG)
    logging.getLogger('PostProcessRemoveBond').setLevel(logging.DEBUG)
    logging.getLogger('PostProcessJoinParticles').setLevel(logging.DEBUG)

    sock = espressopp.tools.vmd.connect(system)

# Reset the integrator step.
    integrator.step = 0
# Run system with non-capped potentials, no thermostat, fixed LJ epsilon
    time0 = time.time()
    traj_file.dump(0, 0)
    print('Running serious simulation %s %s in T=%s' % (
        args.loops*args.steps, 'steps', conf.T))
    system_analysis.info()
    for k in range(args.loops):
        integrator.run(1)
        system_analysis.info()
        traj_file.dump(k*args.steps, k*args.steps*conf.dt)
        dump_topol.update()
        traj_file.flush()
        espressopp.tools.vmd.imd_positions(system, sock)
    dump_topol.update()

    traj_file.close()
    espressopp.tools.analyse.final_info(system, integrator, verletList, time0, time.time())
    a_a_bonds = fpl_a_a.getBonds()
    a_a_a_angles = tpl_a_a.getTriples()
    import cPickle
    cPickle.dump([a_a_bonds, a_a_a_angles], open('topol.pck', 'wb'))
    print('Created {} bonds'.format(len(a_a_bonds[0])))
    print('Created {} angles'.format(len(a_a_a_angles[0])))

if __name__ == '__main__':
    main()
