#  Copyright (C) 2016
#      Jakub Krajniak (jkrajniak at gmail.com)
#
#  This file is part of ESPResSo++.
#
#  ESPResSo++ is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  ESPResSo++ is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program.  If not, see <http://www.gnu.org/licenses/>.

import espressopp
import ConfigParser
import re


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
    reactant_list['type_1']['new_type'] = products[0]['name']
    reactant_list['type_2']['new_type'] = products[1]['name']

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
    reactant_list['type_1']['new_type'] = products[0]['name']
    reactant_list['type_2']['new_type'] = products[1]['name']

    return reactant_list


def process_reaction(reaction):
    reaction = dict(reaction)

    group = reaction['group']
    data = {
        'rate': float(reaction['rate']),
        'cutoff': float(reaction['cutoff']),
        'intramolecular': eval(reaction.get('intramolecular', 'False')),
        'intraresidual': eval(reaction.get('intraresidual', 'False'))
        }

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
    return {
        'interval': int(cfg['interval']),
        'nearest': bool(cfg.get('nearest', False)),
        'bond_limit': int(cfg.get('bond_limit', -1))
    }


def process_group(cfg):
    cfg = dict(cfg)
    return_gr = {'potential': cfg['potential'],
                 'potential_options': dict(
                     [s.split('=') for s in cfg['potential_options'].split(',')]),
                 'reaction_list': [],
                 'extensions': {}
                 }
    if 'extensions' in cfg:
        return_gr['extensions'] = {s.strip(): None for s in cfg['extensions'].split(',')}
    return return_gr


def process_extension(cfg):
    cfg = dict(cfg)
    ret = {'class': cfg['ext_type']}
    del cfg['ext_type']
    ret['options'] = cfg
    return ret


def parse_config(input_file):
    parser = ConfigParser.SafeConfigParser()
    parser.read(input_file)

    config = {'general': None, 'reactions': {}}
    extensions = {}

    for s in parser.sections():
        if s == 'general':
            config['general'] = process_general(parser.items(s))
        elif s.startswith('group_'):
            group_name = s.replace('group_', '').strip()
            group_opt = None
            if group_name not in config['reactions']:
                group_opt = process_group(parser.items(s))
                config['reactions'][group_name] = group_opt
                for ext in group_opt['extensions']:
                    group_opt['extensions'][ext] = extensions[ext]
        elif s.startswith('ext_'):
            name = s.replace('ext_', '').strip()
            properties = process_extension(parser.items(s))
            if name in extensions:
                raise RuntimeError('Name of extension already exists')
            extensions[name] = properties
        elif s.startswith('reaction_'):
            group_name, data = process_reaction(parser.items(s))
            if group_name not in config['reactions']:
                raise RuntimeError(
                    'Wrong order, first reaction groups and then referring reactions')
            config['reactions'][group_name]['reaction_list'].append(data)
    return config


class SetupReactions:
    def __init__(self, system, vl, topol, topol_manager, config):
        self.system = system
        self.vl = vl
        self.topol = topol
        self.tm = topol_manager
        self.cfg = config
        self.name2type = topol.atomsym_atomtype

    def setup_reaction(self, chem_reaction, fpl):
        rl = chem_reaction['reactant_list']
        if not chem_reaction['active']:
            return None
        if chem_reaction['reverse']:
            r_class = espressopp.integrator.DissociationReaction
        else:
            r_class = espressopp.integrator.Reaction
        rt1 = rl['type_1']['name']
        rt2 = rl['type_2']['name']
        r = r_class(
            type_1=self.name2type[rl['type_1']['name']],
            type_2=self.name2type[rl['type_2']['name']],
            delta_1=int(rl['type_1']['delta']),
            delta_2=int(rl['type_2']['delta']),
            min_state_1=int(rl['type_1']['min']),
            max_state_1=int(rl['type_1']['max']),
            min_state_2=int(rl['type_2']['min']),
            max_state_2=int(rl['type_2']['max']),
            rate=float(chem_reaction['rate']),
            fpl=fpl,
            cutoff=float(chem_reaction['cutoff'])
        )
        print('Setup reaction: {}({})-{}({})'.format(
            rt1, self.name2type[rt1], rt2, self.name2type[rt2]))
        if not chem_reaction['reverse']:
            r.intramolecular = bool(chem_reaction['intramolecular'])
            r.intraresidual = bool(chem_reaction['intraresidual'])
        if 'min_cutoff' in chem_reaction:
            r.min_cutoff = float(chem_reaction['min_cutoff'])
        if 'diss_rate' in chem_reaction:
            r.diss_rate = float(chem_reaction['diss_rate'])
        if 'active' in chem_reaction:
            r.active = chem_reaction['active']

        # Change type if necessary.
        if (rl['type_1']['name'] != rl['type_1']['new_type'] or
                rl['type_2']['name'] != rl['type_2']['new_type']):
            r_pp = espressopp.integrator.PostProcessChangeProperty()
            t1_old = self.name2type[rl['type_1']['name']]
            t1_new = self.name2type[rl['type_1']['new_type']]
            if t1_old != t1_new:
                print('Reaction: {}-{}, change type {}->{}'.format(rt1, rt2, t1_old, t1_new))
                new_property = self.topol.gt.atomtypes[rl['type_1']['new_type']]
                r_pp.add_change_property(
                    t1_old,
                    espressopp.ParticleProperties(
                        t1_new, new_property['mass'],
                        new_property['charge']))
            t2_old = self.name2type[rl['type_2']['name']]
            t2_new = self.name2type[rl['type_2']['new_type']]
            if t2_old != t2_new:
                print('Reaction: {}-{}, change type {}->{}'.format(rt1, rt2, t2_old, t2_new))
                new_property = self.topol.gt.atomtypes[rl['type_2']['new_type']]
                r_pp.add_change_property(
                    t2_old,
                    espressopp.ParticleProperties(
                        t2_new, new_property['mass'],
                        new_property['charge']))
            r.add_postprocess(r_pp)

        return r

    def _prepare_group_postprocess(self, cfg):
        pps = []

        def _cfg_post_process_change_neighbour(cfg):
            pp = espressopp.integrator.PostProcessChangeNeighboursProperty(self.tm)
            type_transfers = [
                x.split('->') for x in cfg['type_transfers'].split(',')]
            for old_type, new_type in type_transfers:
                old_type, nb_level = old_type.split(':')
                nb_level = int(nb_level)
                if old_type != new_type:
                    print('Change property {}->{} nb={} and {}'.format(
                        old_type, new_type, nb_level, nb_level+1))
                    t1_old = self.name2type[old_type]
                    t1_new = self.name2type[new_type]
                    new_property = self.topol.gt.atomtypes[new_type]
                    pp.add_change_property(
                        t1_old,
                        espressopp.ParticleProperties(
                            t1_new, new_property['mass'], new_property['charge']),
                        nb_level
                    )
                    pp.add_change_property(
                        t1_old,
                        espressopp.ParticleProperties(
                            t1_new, new_property['mass'], new_property['charge']),
                        nb_level+1
                    )
            return pp

        class_to_cfg = {
            'ChangeNeighboursProperty': _cfg_post_process_change_neighbour
        }
        for pp_cfg in cfg.values():
            cfg_setup = class_to_cfg.get(pp_cfg['class'])
            pps.append(cfg_setup(pp_cfg['options']))

        return pps

    def setup_reactions(self):
        """Setup reactions."""
        self.ar_interval = int(self.cfg['general']['interval'])
        ar = espressopp.integrator.ChemicalReaction(
            self.system,
            self.vl,
            self.system.storage,
            self.tm,
            self.ar_interval)
        ar.nearest_mode = self.cfg['general']['nearest']
        ar.bond_limit = self.cfg['general']['bond_limit']

        fpls = []

        for group_name, reaction_group in self.cfg['reactions'].items():
            print('Setting reaction group {}'.format(group_name))

            # Setting the interaction for the pairs created by this reaction group.
            fpl = espressopp.FixedPairList(self.system.storage)
            fpls.append(fpl)
            pot_class = eval('espressopp.interaction.{}'.format(reaction_group['potential']))
            # Convert if it's possible, values for float
            pot_options = {}
            for k, v in reaction_group['potential_options'].items():
                try:
                    pot_options[k] = float(v)
                except ValueError:
                    pot_options[k] = v
            print('Setting potential for bond with class {}, options {}'.format(
                reaction_group['potential'], reaction_group['potential_options']))
            potential = pot_class(**pot_options)
            interaction = eval('espressopp.interaction.FixedPairList{}'.format(
                reaction_group['potential']))(self.system, fpl, potential)
            self.system.addInteraction(interaction, 'fpl_{}'.format(group_name))

            # Setting the post process extensions.
            extensions = self._prepare_group_postprocess(reaction_group['extensions'])

            print('Setting chemical reactions in group')
            for chem_reaction in reaction_group['reaction_list']:
                r = self.setup_reaction(chem_reaction, fpl)
                if r is not None:
                    for pp in extensions:
                        r.add_postprocess(pp)
                    ar.add_reaction(r)
        return ar, fpls