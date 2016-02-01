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
    re_reactant = re.compile(r'(?P<name>\w+)\((?P<min>\d+),\s+(?P<max>\d+)\)')
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


def process_reaction(reaction):
    reaction = dict(reaction)
    potential_options = dict([s.split('=') for s in reaction['potential_options'].split(',')])

    return (reaction['group'], {
        'rate': reaction['rate'],
        'cutoff': reaction['cutoff'],
        'potential': 'espressopp.interaction.{}'.format(reaction['potential']),
        'potential_options': potential_options,
        'reactant_list': parse_equation(reaction['reaction'])
        })


def process_general(cfg):
    cfg = dict(cfg)
    return {'interval': int(cfg['interval'])}


def parse_config(input_file):
    parser = ConfigParser.SafeConfigParser()
    parser.read(input_file)

    config = {'general': None, 'reactions': {}}

    for s in parser.sections():
        if s == 'general':
            config['general'] = process_general(parser.items(s))
        else:
            group_name, data = process_reaction(parser.items(s))
            if group_name not in config['reactions']:
                config['reactions'][group_name] = []
            config['reactions'][group_name].append(data)

    return config


def setup_reactions(system, verletlist, config):

    ar = espressopp.integrator.ChemicalReaction(
        system, verletlist, system.storage, config['general']['interval'])

    for reaction_group in config['reactions']:
        print reaction_group

    return ar
