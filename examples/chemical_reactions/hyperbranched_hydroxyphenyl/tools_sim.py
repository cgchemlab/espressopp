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

import ast
import collections
import os

import espressopp  # noqa
import random

import tools as general_tools

__doc__ = 'The tools for the simulation.'


Molecule = collections.namedtuple('Molecule', ['pid', 'pos', 'mass', 'type'])


def warmup(system,
           integrator,
           number,
           potential_pairs=None,
           total_force_capping=1000000.0):
    """Warmup the system.

    Args:
        system: The system object.
        integrator: The integrator object.
        number: The number of steps of the warm up phase.
        potential_pairs: The list of tuples with the particle types potentials.
        total_force_capping: The maximum force capping.
    """
    if potential_pairs is None:
        potential_pairs = [(0, 0)]

    print('\nStarting warmup')

    org_dt = integrator.dt
    integrator.dt = pow(10, -4)

    N = 100
    warmup_steps = number / N
    d_force = total_force_capping / warmup_steps

    print('warmup with dt = %e' % integrator.dt)
    print('max_force = %d' % total_force_capping)
    print('number = %d' % number)
    print('warmup time = %e s' % (warmup_steps*integrator.dt))
    print('d_force = %e' % d_force)

    potentials = {}
    interaction = system.getInteractionByLabel('lj')
    for pair in potential_pairs:
        pot = interaction.getPotential(*pair)
        potentials[pair] = [
            pot,
            pot.sigma,
            pot.epsilon,
            pot.epsilon / warmup_steps,
            pot.sigma / warmup_steps
            ]

    force_capping = espressopp.integrator.CapForce(system, d_force)
    integrator.addExtension(force_capping)

    for k in range(1, warmup_steps+1):
        # Update the sigma and the epsilon.
        for pair, pot in potentials.iteritems():
            pot[0].epsilon = k*pot[3]
            pot[0].sigma = k*pot[4]
            interaction.setPotential(pair[0], pair[1], pot[0])
        # e_pots.append(EPot.compute())
        # temps.append(T.compute())
        # print k, numpy.std(e_pots), numpy.std(temps)
        integrator.run(N)
        force_capping.setAbsCapForce(k*d_force)

    # Restore parameters.
    for pair, pot in potentials.iteritems():
        pot[0].sigma = pot[1]
        pot[0].epsilon = pot[2]
        interaction.setPotential(pair[0], pair[1], pot[0])

    force_capping.disconnect()

    print('Force capping disconnected')

    for k in range(2*warmup_steps):
        integrator.run(N)

    integrator.dt = org_dt
    integrator.step = 0
    print("warmup finished")


def warmup_capped(system, integrator, verletList, rc, potential_matrix,
                  sigma_start, warmup_loops, warmup_steps):
    LJ = espressopp.interaction.VerletListLennardJones(verletList)
    for type1, type2, sigma, epsilon in potential_matrix:
        LJ.setPotential(
            type1=type1,
            type2=type2,
            potential=espressopp.interaction.LennardJones(
                sigma=sigma,
                epsilon=epsilon,
                cutoff=rc
            ))
    system.addInteraction(LJ, 'lj_warmup')

    old_dt = integrator.dt
    integrator.dt = 0.0001

    espressopp.tools.analyse.info(system, integrator, kb=system.kb)

    print('Running with capped potential, increasing sigma, loops {}, steps {}'.format(
        warmup_loops, warmup_steps))
    # Run system with capped potentials, thermostat and increasing LJ epsilon
    for k in range(warmup_loops):
        for type1, type2, sigma, epsilon in potential_matrix:
            LJ.setPotential(
                type1=type1,
                type2=type2,
                potential=espressopp.interaction.LennardJones(
                    sigma=sigma_start + (sigma-sigma_start)*k*1.0/(warmup_loops-1),
                    epsilon=epsilon,
                    cutoff=rc
                ))
        integrator.run(warmup_steps)
        if integrator.step % 100 == 0:
            espressopp.tools.analyse.info(system, integrator, kb=system.kb)

    # Remove LJ Capped potential
    system.removeInteractionByLabel('lj_warmup')
    integrator.dt = old_dt


def combination(sig_1, eps_1, sig_2, eps_2, cr):
    if cr == 2:
        sig = 0.5*(sig_1 + sig_2)
        eps = (eps_1*eps_2)**(1.0/2.0)
    else:
        sig = (sig_1*sig_2)**(1.0/2.0)
        eps = (eps_1*eps_2)**(1.0/2.0)

    return sig, eps


def setNonbondedInteractions(system, gt, vl, lj_cutoff, tab_cutoff=None):  #NOQA
    defaults = gt.gt.defaults
    atomparams = gt.gt.atomtypes
    atomsym_atomtype = gt.atomsym_atomtype

    if tab_cutoff is None:
        tab_cutoff = lj_cutoff

    combinationrule = int(defaults['combinationrule'])
    print('Settings up LJ interactions')
    type_pairs = set()
    for type_1 in atomsym_atomtype:
        for type_2 in atomsym_atomtype:
            type_pairs.add(tuple(sorted([type_1, type_2])))

    lj_interaction = espressopp.interaction.VerletListLennardJones(vl)
    tab_interaction = espressopp.interaction.VerletListTabulated(vl)

    # Special case for MultiTabulated
    cr_multi = collections.defaultdict(list)
    cr_observs = {}

    print('Number of pairs: {}'.format(len(type_pairs)))
    for type_1, type_2 in type_pairs:
        t1 = atomsym_atomtype[type_1]
        t2 = atomsym_atomtype[type_2]
        param = gt.gt.nonbond_params.get((type_1, type_2))
        table_name = None
        cr_type = None
        cr_min, cr_max = 0, 0
        cr_total = 0
        sig_1, eps_1, sig_2, eps_2 = 0, 0, 0, 0
        sig, eps = -1, -1
        if param:
            print('Using defined non-bonded cross params')
            func = param['func']
            if func == 1:
                sig = float(param['params'][0])
                eps = float(param['params'][1])
            elif func == 8:
                table_name = 'table_{}_{}.xvg'.format(type_1, type_2)
            elif func == 1:
                sig_1, eps_1 = atomparams[type_1]['sig'], atomparams[type_1]['eps']
                sig_2, eps_2 = atomparams[type_2]['sig'], atomparams[type_2]['eps']
                sig, eps = combination(sig_1, eps_1, sig_2, eps_2, combinationrule)
            elif func == 9:
                tab_name = 'table_{}_{}.xvg'.format(param['params'][1], param['params'][0])
                cr_type = atomsym_atomtype[param['params'][2]]
                cr_total = int(param['params'][3])
                cr_min = float(param['params'][4])
                cr_max = float(param['params'][5])
                cr_default = bool(int(param['params'][6])) if len(param['params']) > 6 else False
                if (cr_type, cr_total) not in cr_observs:
                    cr_observs[(cr_type, cr_total)] = espressopp.analysis.ChemicalConversion(
                        system, cr_type, cr_total)
                cr_multi[(t1, t2)].append([
                    cr_observs[(cr_type, cr_total)],
                    tab_name,
                    cr_min,
                    cr_max,
                    cr_default])
        else:
            sig_1, eps_1 = atomparams[type_1]['sig'], atomparams[type_1]['eps']
            sig_2, eps_2 = atomparams[type_2]['sig'], atomparams[type_2]['eps']
            sig, eps = combination(sig_1, eps_1, sig_2, eps_2, combinationrule)
        # Standard interaction.
        if sig > 0 and eps > 0:
            print('Set lj potential {}-{}, eps={}, sig={}'.format(type_1, type_2, eps, sig))
            ljpot = espressopp.interaction.LennardJones(
                epsilon=eps, sigma=sig, cutoff=lj_cutoff)
            lj_interaction.setPotential(type1=t1, type2=t2, potential=ljpot)
        elif table_name is not None:
            print('Set tab potential {}-{}: {}'.format(type_1, type_2, table_name))
            espp_tab_name = '{}.pot'.format(table_name.replace('.xvg', ''))
            if not os.path.exists(espp_tab_name):
                print('Convert {} to {}'.format(table_name, espp_tab_name))
                espressopp.tools.convert.gromacs.convertTable(table_name, espp_tab_name)
            tab_interaction.setPotential(
                type1=t1, type2=t2,
                potential=espressopp.interaction.Tabulated(
                    itype=2, filename=espp_tab_name, cutoff=tab_cutoff))

    system.addInteraction(lj_interaction, 'lj')
    system.addInteraction(tab_interaction, 'lj-tab')

    # Mixed tabulated potentials.
    if cr_multi:
        multi_tab_interaction = espressopp.interaction.VerletListMultiTabulated(vl)
        for (mt1, mt2), data in cr_multi.items():
            mp_tab = espressopp.interaction.MultiTabulated(cutoff=tab_cutoff)
            for cr_obs, tab_name, cr_min, cr_max, cr_default in data:
                espp_tab_name = '{}.pot'.format(tab_name.replace('.xvg', ''))
                if not os.path.exists(espp_tab_name):
                    print('Convert {} to {}'.format(tab_name, espp_tab_name))
                    espressopp.tools.convert.gromacs.convertTable(tab_name, espp_tab_name)
                mp_tab.register_table(espp_tab_name, 2, cr_obs, cr_min, cr_max, cr_default)
            print('Set multi tabulated potential {}-{}'.format(mt1, mt2))
            multi_tab_interaction.setPotential(
                type1=mt1, type2=mt2, potential=mp_tab)
    system.addInteraction(multi_tab_interaction, 'lj-mtab')
    return cr_observs


def setBondInteractions(system, gt):
    fpl = espressopp.FixedPairList(system.storage)
    fpl.addBonds(gt.bonds)
    tab_interaction = espressopp.interaction.FixedPairListTypesTabulated(system, fpl)
    for (t1, t2), param in gt.bondparams.items():
        if param['func'] != 8:
            raise RuntimeError('Wrong func type, only tabulated supported')
        espp_tab_name = 'table_b{}.pot'.format(param['params'][0])
        tab_name = 'table_b{}.xvg'.format(param['params'][0])
        if not os.path.exists(espp_tab_name):
            print('Convert {} to {}'.format(tab_name, espp_tab_name))
            espressopp.tools.convert.gromacs.convertTable(tab_name, espp_tab_name)
        tab_interaction.setPotential(
            type1=t1, type2=t2,
            potential=espressopp.interaction.Tabulated(2, espp_tab_name))
    system.addInteraction(tab_interaction, 'bonds')
    return fpl, tab_interaction


def setAngleInteractions(system, gt):
    fpl = espressopp.FixedTripleList(system.storage)
    fpl.addTriples(gt.angles)
    tab_interaction = espressopp.interaction.FixedTripleListTypesTabulatedAngular(system, fpl)
    for (t1, t2, t3), param in gt.angleparams.items():
        if param['func'] != 8:
            raise RuntimeError('Wrong func type, only tabulated supported')
        espp_tab_name = 'table_a{}.pot'.format(param['params'][0])
        tab_name = 'table_a{}.xvg'.format(param['params'][0])
        if not os.path.exists(espp_tab_name):
            print('Convert {} to {}'.format(tab_name, espp_tab_name))
            espressopp.tools.convert.gromacs.convertTable(tab_name, espp_tab_name)
        tab_interaction.setPotential(
            type1=t1, type2=t2, type3=t3,
            potential=espressopp.interaction.TabulatedAngular(2, espp_tab_name))
    system.addInteraction(tab_interaction, 'angles')
    return fpl, tab_interaction


def setDihedralInteractions(system, gt):
    fpl = espressopp.FixedQuadrupleList(system.storage)
    fpl.addQuadruples(gt.dihedrals)
    tab_interaction = espressopp.interaction.FixedQuadrupleListTypesTabulatedDihedral(system, fpl)
    for (t1, t2, t3, t4), param in gt.dihedralparams.items():
        if param['func'] != 8:
            raise RuntimeError('Wrong func type, only tabulated supported')
        espp_tab_name = 'table_d{}.pot'.format(param['params'][0])
        tab_name = 'table_d{}.xvg'.format(param['params'][0])
        if not os.path.exists(espp_tab_name):
            print('Convert {} to {}'.format(tab_name, espp_tab_name))
            espressopp.tools.convert.gromacs.convertTable(tab_name, espp_tab_name)
        tab_interaction.setPotential(
            type1=t1, type2=t2, type3=t3, type4=t4,
            potential=espressopp.interaction.TabulatedDihedral(2, espp_tab_name))
    system.addInteraction(tab_interaction, 'dihedrals')
    return fpl, tab_interaction


def genParticleList(coordinate, topol):
    """Generates particle list
    Args:
        coordinate: The input coordinate file.
        topol: Topology file.
    Returns:
        List of property names and particle list.
    """
    props = ['id', 'type', 'pos', 'mass', 'q', 'res_id', 'state']
    particle_list = []

    for atom_id in sorted(coordinate.atoms):
        data = coordinate.atoms[atom_id]
        top_data = topol.atomparams['{}-{}'.format(data.chain_name, data.name)]
        particle_list.append(
            [atom_id,
             top_data['type_id'],
             espressopp.Real3D(data.position),
             top_data['mass'],
             top_data['charge'],
             data.chain_idx,
             top_data.get('state', 0)]
        )

    return props, particle_list


def setTopologyManager(input_conf, tm, bondedint, angleint, dihedralint, pairint):
    pid2type = {p+1: input_conf.types[p] for p in xrange(len(input_conf.types))}
    # Set input tuples. That does not change but serves as an input.
    for bi in bondedint.values():
        tm.observe_tuple(bi.getFixedPairList())

    # Register angle, dihedral, pair fixed pair lists.
    for ai in angleint.values():
        fl = ai.getFixedTripleList()
        triples = fl.getTriples()
        type_triples = {tuple(map(pid2type.get, m)) for t in triples for m in t}
        for types in type_triples:
            tm.register_triplet(fl, *types)

    # Register dihedrals
    for di in dihedralint.values():
        fl = di.getFixedQuadrupleList()
        quadruples = fl.getQuadruples()
        type_quadruples = {tuple(map(pid2type.get, m)) for t in quadruples for m in t}
        for types in type_quadruples:
            tm.register_quadruplet(fl, *types)

    for pi in pairint.values():
        fl = pi.getFixedPairList()
        pairs = fl.getPairs()
        type_pairs = {tuple(map(pid2type.get, m)) for t in pairs for m in t}
        for types in type_pairs:
            tm.register_tuple(fl, *types)


def setDynamicExcludeList(dyn, bondedint, angleint, dihedralint, pairint):
    for bi in bondedint.values():
        dyn.observe_tuple(bi.getFixedPairList())
    for ai in angleint.values():
        dyn.observe_triple(ai.getFixedTripleList())
    for di in dihedralint.values():
        dyn.observe_quadruple(di.getFixedQuadrupleList())
    for pi in pairint.values():
        dyn.observe_tuple(pi.getFixedPairList())


def _args():
    parser = general_tools.MyArgParser(description='Runs classical MD simulation',
                                       fromfile_prefix_chars='@')
    parser.add_argument('--conf', required=True, help='Input .gro coordinate file')
    parser.add_argument('--top', '--topology', required=True, help='Topology file',
                        dest='top')
    parser.add_argument('--node_grid')
    parser.add_argument('--skin', type=float, default=0.16,
                        help='Skin value for Verlet list')
    parser.add_argument('--coord', help='Input coordinate h5md file')
    parser.add_argument('--coord_frame', default=-1, type=int,
                        help='Time frame of input coordinate h5md file')
    parser.add_argument('--run', type=int, default=10000,
                        help='Number of simulation steps')
    parser.add_argument('--int_step', default=1000, type=int, help='Steps in integrator')
    parser.add_argument('--rng_seed', type=int, help='Seed for RNG', required=False,
                        default=random.randint(1000, 10000))
    parser.add_argument('--output_prefix',
                        default='sim', type=str,
                        help='Prefix for output files')
    parser.add_argument('--output_file',
                        default='trjout.h5', type=str,
                        help='Name of output trajectory file')
    parser.add_argument('--thermostat',
                        default='lv',
                        choices=('lv', 'vr'),
                        help='Thermostat to use, lv: Langevine, vr: Stochastic velocity rescale')
    parser.add_argument('--barostat', default='lv', choices=('lv', 'br'),
                        help='Barostat to use, lv: Langevine, br: Berendsen')
    parser.add_argument('--barostat_tau', default=5.0, type=float,
                        help='Tau parameter for Berendsen barostat')
    parser.add_argument('--barostat_mass', default=50.0, type=float,
                        help='Mass parameter for Langevin barostat')
    parser.add_argument('--barostat_gammaP', default=1.0, type=float,
                        help='gammaP parameter for Langevin barostat')
    parser.add_argument('--thermostat_gamma', type=float, default=5.0,
                        help='Thermostat coupling constant')
    parser.add_argument('--temperature', default=458.0, type=float, help='Temperature')
    parser.add_argument('--pressure', help='Pressure', type=float)
    parser.add_argument('--trj_collect', default=1000, type=int,
                        help='Collect trajectory every (step)')
    parser.add_argument('--energy_collect', default=1000, type=int,
                        help='Collect energy every (step)')
    parser.add_argument('--dt', default=0.001, type=float,
                        help='Integrator time step')
    parser.add_argument('--lj_cutoff', default=1.2, type=float,
                        help='Cutoff of atomistic non-bonded interactions')
    parser.add_argument('--cg_cutoff', default=1.4, type=float,
                        help='Cuoff of coarse-grained non-bonded interactions')
    parser.add_argument('--coulomb_epsilon1', default=1.0, type=float,
                        help='Epsilon_1 for coulomb interactions')
    parser.add_argument('--coulomb_epsilon2', default=80.0, type=float,
                        help='Epsilon_2 for coulomb interactions')
    parser.add_argument('--coulomb_kappa', default=1.0, type=float,
                        help='Kappa paramter for coulomb interactions')
    parser.add_argument('--table_groups', default='A,B',
                        help='Name of CG groups to read from tables')
    parser.add_argument('--initial_step', default=0,
                        help='Initial integrator step (useful for continue simulation',
                        type=int)
    parser.add_argument('--reactions', default='reaction.cfg',
                        help='Configuration file with chemical reactions')
    parser.add_argument('--debug', default=None, help='Turn on logging mechanism')
    parser.add_argument('--start_ar', default=0, type=int, help='When to start chemical reactions')
    parser.add_argument('--interactive', default=False, type=ast.literal_eval,
                        help='Run interactive mode')
    parser.add_argument('--store_species', default=False, type=ast.literal_eval,
                        help='Store particle types')
    parser.add_argument('--store_state', default=True, type=ast.literal_eval,
                        help='Store chemical state')
    parser.add_argument('--store_lambda', default=False, type=ast.literal_eval,
                        help='Store lambda parameter')

    return parser
