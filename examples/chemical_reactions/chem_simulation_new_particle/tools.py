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


import argparse
import files_io
import sys

__doc__ = "Tool functions."


class MyArgParser(argparse.ArgumentParser):
    def __init__(self, *args, **kwargs):
        super(MyArgParser, self).__init__(*args, **kwargs)

    def convert_arg_line_to_args(self, line):
        for arg in line.split():
            t = arg.strip()
            if not t:
                continue
            if t.startswith('#'):
                break
            if not t.startswith('--'):
                t = '--{}'.format(t)
            yield t

    @staticmethod
    def save_to_file(output_file, namespace):
        """Saves arguments to file so it can be read again.

        Args:
            output_file: The string with the name of output file.
            namespace: The namespace with arguements.
        """
        with open(output_file, "w") as of:
            keys = sorted(namespace.__dict__.keys())
            for k in keys:
                v = namespace.__dict__[k]
                if v is not None:
                    of.write('{}={}\n'.format(k, v))


def get_graph(settings):
    """Build graph based on settings file. Useful for GROMACS."""
    gro = files_io.GROFile(settings.cg_configuration['file'])
    gro.read()

    g = nx.Graph(box=gro.box)
    for at_id, at in gro.atoms.iteritems():
        g.add_node(
            at_id,
            name=at.name,
            res_id=at.chain_idx,
            position=at.position,
            chain_name=at.chain_name)

    # Adding edges
    for mol in gro.chains:
        try:
            cg_bonds = settings.cg_molecules[mol].molecule_topology.get('bond')
        except KeyError:
            print(('\nError:\nMolecule \'{}\' not found in input CG trajectory'
                   '(valid molecule\' names: {})\nExit, nothing to do.'
                   ).format(mol, settings.cg_molecules.keys()))
            sys.exit(1)
        if cg_bonds:
            for chain_idx in gro.chains[mol]:
                for bond_name, bond_def in cg_bonds.iteritems():
                    for b1, b2 in bond_def['list']:
                        a1 = gro.chains[mol][chain_idx][b1]
                        a2 = gro.chains[mol][chain_idx][b2]
                        g.add_edge(a1.atom_id, a2.atom_id, params=bond_def['params'])
    # Update degree
    for n_id in g.node:
        g.node[n_id]['degree'] = g.degree(n_id)

    return g
