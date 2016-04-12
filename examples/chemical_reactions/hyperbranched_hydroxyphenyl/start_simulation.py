#!/usr/bin/env python
"""
Copyright (C) 2015-2016 Jakub Krajniak <jkrajniak@gmail.com>

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

import espressopp  # NOQA
import math  # NOQA
try:
    import MPI
except ImportError:
    from mpi4py import MPI
import time
import logging
import random
import shutil

import gromacs_topology_new
import files_io
import reaction_parser
import tools_sim

# GROMACS units, kJ/mol K
kb = 0.0083144621

h5md_group = 'atoms'

__doc__ = 'Run GROMACS-like simulation'

# Do not to modify lines below.


def sort_trajectory(trj, ids):
    """Performs sorting on HDF5 file. It is required because by default, H5MD file
    can be in unsorted state and only /particles/{}/id/value inform about particle. id.

    Args:
        trj: The input HDF5 Dataset to sort.
        ids: The input ids Dataset with particle ids for every time step.

    Returns:
        Sorted numpy array.
    """
    print('Sorting trajectory')
    idd = [
        x[1] for x in sorted([(p_id, col_id) for col_id, p_id in enumerate(ids)],
                             key=lambda y: (True, y[0]) if y[0] == -1 else (False, y[0]))
    ]
    return trj[idd]


def main():  #NOQA
    args = tools_sim._args().parse_args()

    tools_sim._args().save_to_file('{}params.out'.format(args.output_prefix), args)

    if args.debug:
        for s in args.debug.split(','):
            print('Activating logger {}'.format(s))
            logging.getLogger(s.strip()).setLevel(logging.DEBUG)

    lj_cutoff = args.lj_cutoff
    cg_cutoff = args.cg_cutoff
    max_cutoff = max([lj_cutoff, cg_cutoff])
    dt = args.dt

    time0 = time.time()

    gt = gromacs_topology_new.GromacsTopology(args.top)
    gt.read()

    input_conf = files_io.GROFile(args.conf)
    input_conf.read()

    box = input_conf.box
    print('Setup simulation...')

    # Tune simulation parameter according to arguments
    integrator_step = args.int_step
    sim_step = args.run / integrator_step

    if args.skin:
        skin = args.skin

    rng_seed = args.rng_seed
    if not args.rng_seed:
        rng_seed = random.randint(10, 1000000)

    print('Skin: {}'.format(skin))
    print('RNG Seed: {}'.format(rng_seed))

    part_prop, particle_list = tools_sim.genParticleList(input_conf, gt)
    print('Reads {} particles with properties {}'.format(len(particle_list), part_prop))

    density = sum(x[3] for x in particle_list)*1.6605402 / (box[0] * box[1] * box[2])
    print('Density: {} kg/m^3'.format(density))
    print('Box: {} nm'.format(box))

    print gt.atomsym_atomtype

    # Generate velocity.
    print('Generating velocities from Maxwell-Boltzmann distribution T={}'.format(
        args.temperature))
    vx, vy, vz = espressopp.tools.velocities.gaussian(
        args.temperature, len(particle_list), [x[3]*1.6605402 for x in particle_list],
        kb=kb)
    part_prop.append('v')
    for i, p in enumerate(particle_list):
        p.append(espressopp.Real3D(vx[i], vy[i], vz[i]))

    system = espressopp.System()
    system.rng = espressopp.esutil.RNG(rng_seed)
    system.bc = espressopp.bc.OrthorhombicBC(system.rng, box)
    system.skin = skin
    if args.node_grid:
        nodeGrid = map(int, args.node_grid.split(','))
    else:
        nodeGrid = espressopp.tools.decomp.nodeGrid(MPI.COMM_WORLD.size)
    print('Number of nodes {}, node-grid: {}'.format(
        MPI.COMM_WORLD.size, nodeGrid))
    cellGrid = espressopp.tools.decomp.cellGrid(box, nodeGrid, max_cutoff, skin)

    print('Cell grid: {}'.format(cellGrid))

    system.storage = espressopp.storage.DomainDecomposition(system, nodeGrid, cellGrid)
    integrator = espressopp.integrator.VelocityVerlet(system)
    integrator.dt = dt
    system.integrator = integrator

    system.storage.addParticles(particle_list, *part_prop)
    system.storage.decompose()

# In the case of butane is very easy to do
    dynamic_exclusion_list = espressopp.DynamicExcludeList(integrator, gt.exclusions)
    print('Excluded pairs from LJ interaction: {}'.format(len(gt.exclusions)))

# Exclude all bonded interaction from the lennard jones
    verletlist = espressopp.VerletList(
        system,
        cutoff=max_cutoff,
        exclusionlist=dynamic_exclusion_list
        )

# define the potential, interaction_id = 0
    print('Bonds: {}'.format(len(gt.bonds)))
    print('Angles: {}'.format(len(gt.angles)))
    print('Dihedrals: {}'.format(len(gt.dihedrals)))

# Define the thermostat
    temperature = args.temperature*kb
    print('Temperature: {}, gamma: {}'.format(temperature, args.thermostat_gamma))
    print('Thermostat: {}'.format(args.thermostat))
    if args.thermostat == 'lv':
        thermostat = espressopp.integrator.LangevinThermostat(system)
        thermostat.temperature = temperature
        thermostat.gamma = args.thermostat_gamma
    elif args.thermostat == 'vr':
        thermostat = espressopp.integrator.StochasticVelocityRescaling(system)
        thermostat.temperature = temperature
        thermostat.coupling = args.thermostat_gamma
    else:
        raise Exception('Wrong thermostat keyword: `{}`'.format(args.thermostat))
    integrator.addExtension(thermostat)

    pressure_comp = espressopp.analysis.Pressure(system)
    if args.pressure:
        pressure = args.pressure * 0.060221374  # convert from bars to gromacs units kj/mol/nm^3
        if args.barostat == 'lv':
            print('Barostat: Langevin with P={}, gamma={}, mass={}'.format(
                pressure, 0.5, pow(10, 4)))
            barostat = espressopp.integrator.LangevinBarostat(system, system.rng, temperature)
            barostat.gammaP = args.barostat_gammaP
            barostat.mass = args.barostat_mass
            barostat.pressure = pressure
        elif args.barostat == 'br':
            print('Barostat: Berendsen with P={} and tau={}'.format(pressure, 0.5))
            barostat = espressopp.integrator.BerendsenBarostat(system, pressure_comp)
            barostat.tau = args.barostat_tau
            barostat.pressure = pressure
        integrator.addExtension(barostat)

    print("Decomposing now ...")
    system.storage.decompose()

    # Set potentials.
    cr_observs = tools_sim.setNonbondedInteractions(system, gt, verletlist, lj_cutoff, cg_cutoff)
    static_fpl, b_interaction = tools_sim.setBondInteractions(system, gt)
    #static_ftl, _ = tools_sim.setAngleInteractions(system, gt)
    #static_fql, _ = tools_sim.setDihedralInteractions(system, gt)

    print('Set Dynamic Exclusion lists.')
    dynamic_exclusion_list.observe_tuple(static_fpl)
    #dynamic_exclusion_list.observe_triple(static_ftl)
    #dynamic_exclusion_list.observe_quadruple(static_fql)

    print('Set topology manager')
    topology_manager = espressopp.integrator.TopologyManager(system)
    topology_manager.rebuild()
    topology_manager.observe_tuple(static_fpl)
    topology_manager.initialize_topology()
    topology_manager.register_tuple(static_fpl, 0, 0)
    #for t in gt.angleparams:
    #    topology_manager.register_triplet(static_ftl, *t)
    #for t in gt.dihedralparams:
    #    topology_manager.register_quadruplet(static_fql, *t)
    integrator.addExtension(topology_manager)

    # Set chemical reactions
    fpls = []
    cr_interval = 0
    if args.reactions:
        print('Set chemical reactions from: {}'.format(args.reactions))
        reaction_config = reaction_parser.parse_config(args.reactions)
        sc = reaction_parser.SetupReactions(
            system, verletlist, gt, topology_manager, reaction_config)

        ar, fpls = sc.setup_reactions()
        output_reaction_config = '{}_{}_{}'.format(args.output_prefix, rng_seed, args.reactions)
        print('Save copy of reaction config to: {}'.format(output_reaction_config))
        shutil.copyfile(args.reactions, output_reaction_config)
        integrator.addExtension(ar)
        cr_interval = sc.ar_interval

    for f in fpls:
        topology_manager.observe_tuple(f)
        dynamic_exclusion_list.observe_tuple(f)

    energy_file = '{}_energy_{}.csv'.format(args.output_prefix, rng_seed)
    print('Energy saved to: {}'.format(energy_file))
    system_analysis = espressopp.analysis.SystemMonitor(
        system,
        integrator,
        espressopp.analysis.SystemMonitorOutputCSV(energy_file))
    temp_comp = espressopp.analysis.Temperature(system)
    system_analysis.add_observable('T', temp_comp)
    system_analysis.add_observable(
        'Ekin', espressopp.analysis.KineticEnergy(
            system, temp_comp))
    for label, interaction in sorted(system.getAllInteractions().items()):
        print('System analysis: adding {}'.format(label))
        system_analysis.add_observable(
            label, espressopp.analysis.PotentialEnergy(system, interaction))
    for (cr_type, _), obs in cr_observs.items():
        system_analysis.add_observable(
            'cr_{}'.format(cr_type), obs)
    for fidx, f in enumerate(fpls):
        system_analysis.add_observable(
            'count_{}'.format(fidx), espressopp.analysis.NFixedPairListEntries(system, f))
    #system_analysis.add_observable(
    #    'cnt_fpl', espressopp.analysis.NFixedPairListEntries(system, static_fpl))
    #system_analysis.add_observable(
    #    'cnt_ftl', espressopp.analysis.NFixedTripleListEntries(system, static_ftl))
    #system_analysis.add_observable(
    #    'cnt_fql', espressopp.analysis.NFixedQuadrupleListEntries(system, static_fql))

    ext_analysis = espressopp.integrator.ExtAnalyze(system_analysis, cr_interval)
    integrator.addExtension(ext_analysis)
    print('Configured system analysis')

    print('Configure H5MD trajectory writer')
    h5md_output_file = '{}_{}_traj.h5'.format(args.output_prefix, rng_seed)
    traj_file = espressopp.io.DumpH5MD(
        system,
        h5md_output_file,
        group_name='atoms',
        static_box=False,
        author='XXX',
        email='xxx',
        store_species=True,
        store_state=True)
    traj_file.set_parameters({
        'temperature': args.temperature
    })
    print('Set topology writer')
    dump_topol = espressopp.io.DumpTopology(system, integrator, traj_file)
    for i, f in enumerate(fpls):
        dump_topol.observe_tuple(f, 'chem_bonds_{}'.format(i))

    dump_topol.add_static_tuple(static_fpl, 'bonds')
    dump_topol.dump()
    dump_topol.update()
    ext_dump = espressopp.integrator.ExtAnalyze(dump_topol, cr_interval)
    integrator.addExtension(ext_dump)

    print('Reset total velocity')
    total_velocity = espressopp.analysis.TotalVelocity(system)
    total_velocity.reset()

    traj_file.dump(0, 0)

    print('Running {} steps'.format(sim_step*integrator_step))
    system_analysis.dump()
    system_analysis.info()
    for k in range(sim_step):
        integrator.run(integrator_step)
        system_analysis.info()
        total_velocity.reset()
        dump_topol.update()
        traj_file.dump(k*integrator_step, k*integrator_step*dt)
        traj_file.flush()
    else:
        dump_topol.update()
        traj_file.dump(sim_step*integrator_step, sim_step*integrator_step*dt)
        traj_file.close()

    # Saves output file.
    output_gro_file = '{}_{}_confout.gro'.format(args.output_prefix, rng_seed)
    dump_gro = espressopp.io.DumpGRO(
        system, integrator, filename=output_gro_file,
        unfolded=True, append=False)
    dump_gro.dump()
    print('Wrote end configuration to: {}'.format(output_gro_file))

    print('finished!')
    print('total time: {}'.format(time.time()-time0))
    espressopp.tools.analyse.final_info(system, integrator, verletlist, time0, time.time())


if __name__ == '__main__':
    main()
