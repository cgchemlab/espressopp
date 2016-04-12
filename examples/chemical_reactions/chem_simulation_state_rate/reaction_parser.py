#! /usr/bin/env python
#
# Copyright (c) 2016 Jakub Krajniak <jkrajniak@gmail.com>
#
# Distributed under terms of the GNU GPLv3 license.
#

import ConfigParser
import re

import espressopp


def parse_equation(input_string):
    re_reactant = re.compile(r'(?P<name>\w+)\((?P<min>\d+),\s*(?P<max>\d+)\)')
    re_product = re.compile(r'(?P<name>\w+)\((?P<delta>[0-9-]+)\)')

    reactants, products = input_string.split('->')

    reactant_list = {}

    mol_a, mol_b = map(re_reactant.match, map(str.strip, reactants.split('+')))

    reactant_list['type_1'] = mol_a.groupdict()
    reactant_list['type_2'] = mol_b.groupdict()

    products = [re_product.match(x).groupdict() for x in map(str.strip, products.split(':'))]
    reactant_list['type_1']['delta'] = products[0]['delta']
    reactant_list['type_2']['delta'] = products[1]['delta']

    return reactant_list


def parse_reverse_equation(input_string):
    re_reactant = re.compile(r'(?P<name>\w+)\((?P<min>\d+),\s*(?P<max>\d+)\)')
    re_product = re.compile(r'(?P<name>\w+)\((?P<delta>[0-9-]+)\)')

    reactant_list = {}

    reactants, products = map(str.strip, input_string.split('->'))
    mol_a, mol_b = map(re_reactant.match, map(str.strip, reactants.split(':')))

    reactant_list['type_1'] = mol_a.groupdict()
    reactant_list['type_2'] = mol_b.groupdict()
    products = [re_product.match(x).groupdict() for x in map(str.strip, products.split('+'))]
    reactant_list['type_1']['delta'] = products[0]['delta']
    reactant_list['type_2']['delta'] = products[1]['delta']

    return reactant_list


def process_reaction(reaction):
    reaction = dict(reaction)

    group = reaction['group']
    data = {
        'cutoff': float(reaction['cutoff']),
        'intramolecular': eval(reaction.get('intramolecular', 'False')),
        }

    # Handle state dependent reaction rate.
    if 'rates_1' in reaction or 'rates_2' in reaction:
        if 'rate' in reaction:
            raise Exception(
                'Wrong keyword, "rate" cannot be used together with "rates_1" and "rates_2"')
        if 'rates_1' not in reaction or 'rates_2' not in reaction:
            raise Exception(
                'Missing one of keywords, "rates_1" or "rates_2"')
        data['rates_1'] = map(float, reaction['rates_1'].split(','))
        data['rates_2'] = map(float, reaction['rates_2'].split(','))
    elif 'rate' in reaction:
        data['rate'] = float(reaction['rate'])

    try:
        data['reactant_list'] = parse_equation(reaction['reaction'])
        data['reverse'] = False
    except:
        data['reactant_list'] = parse_reverse_equation(reaction['reaction'])
        data['reverse'] = True

    if 'min_cutoff' in reaction:
        data['min_cutoff'] = float(reaction['min_cutoff'])

    if 'diss_rate' in reaction:
        if not data['reverse']:
            raise Exception('Invalid keyword `diss_rate` for non-dissociation reaction')
        data['diss_rate'] = float(reaction['diss_rate'])

    if 'active' in reaction:
        data['active'] = eval(reaction['active'])
    else:
        data['active'] = True

    return (group, data)


def process_general(cfg):
    cfg = dict(cfg)
    return {'interval': int(cfg['interval'])}


def process_group(cfg):
    cfg = dict(cfg)
    return {'potential': cfg['potential'],
            'potential_options': dict([s.split('=') for s in cfg['potential_options'].split(',')]),
            'reaction_list': []}


def parse_config(input_file):
    parser = ConfigParser.SafeConfigParser()
    parser.read(input_file)

    config = {'general': None, 'reactions': {}}

    for s in parser.sections():
        if s == 'general':
            config['general'] = process_general(parser.items(s))
        elif s.startswith('group_'):
            group_name = s.replace('group_', '').strip()
            if group_name not in config['reactions']:
                config['reactions'][group_name] = process_group(parser.items(s))
        else:
            group_name, data = process_reaction(parser.items(s))
            if group_name not in config['reactions']:
                raise RuntimeError(
                    'Wrong order, first reaction groups and then refering reactions')
            config['reactions'][group_name]['reaction_list'].append(data)
    return config


def setup_reactions(system, verletlist, input_conf, config):
    ar = espressopp.integrator.ChemicalReaction(
        system, verletlist, system.storage, int(config['general']['interval']))

    atom_name2atom_type = {
        v['atnum']: k for k, v in input_conf.atomtypeparams.iteritems()
    }
    fpls = []
    rs = []

    for group_name, reaction_group in config['reactions'].items():
        print('Setting reaction group {}'.format(group_name))
        fpl = espressopp.FixedPairList(system.storage)
        fpls.append(fpl)
        pot_class = eval('espressopp.interaction.{}'.format(reaction_group['potential']))
        potential = pot_class(**reaction_group['potential_options'])
        interaction = eval('espressopp.interaction.FixedPairList{}'.format(
            reaction_group['potential']))(system, fpl, potential)
        system.addInteraction(interaction, 'fpl_{}'.format(group_name))
        print('Setting chemical reactions in group')

        for chem_reaction in reaction_group['reaction_list']:
            rl = chem_reaction['reactant_list']
            if chem_reaction['reverse']:
                r_class = espressopp.integrator.DissociationReaction
            else:
                r_class = espressopp.integrator.Reaction
            print('Setting up reaction {}'.format(r_class))

            if 'rates_1' in chem_reaction:
                rates = (chem_reaction['rates_1'], chem_reaction['rates_2'])
            else:
                rates = chem_reaction['rate']

            r = r_class(
                type_1=atom_name2atom_type[rl['type_1']['name']],
                type_2=atom_name2atom_type[rl['type_2']['name']],
                delta_1=int(rl['type_1']['delta']),
                delta_2=int(rl['type_2']['delta']),
                min_state_1=int(rl['type_1']['min']),
                max_state_1=int(rl['type_1']['max']),
                min_state_2=int(rl['type_2']['min']),
                max_state_2=int(rl['type_2']['max']),
                rate=rates,
                fpl=fpl,
                cutoff=float(chem_reaction['cutoff'])
            )
            if not chem_reaction['reverse']:
                r.intramolecular = bool(chem_reaction['intramolecular'])
            if 'min_cutoff' in chem_reaction:
                r.min_cutoff = float(chem_reaction['min_cutoff'])
            if 'diss_rate' in chem_reaction:
                r.diss_rate = float(chem_reaction['diss_rate'])
            if 'active' in chem_reaction:
                r.active = chem_reaction['active']
            ar.add_reaction(r)
    return ar, fpls, rs

if __name__ == '__main__':
    import sys
    print parse_config(sys.argv[1])
    print('Input file is correct')
