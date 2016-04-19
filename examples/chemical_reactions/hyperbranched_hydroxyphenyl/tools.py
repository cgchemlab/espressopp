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
from files_io import TopoAtom

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


def dump_topol(file_name, topol, system, particle_ids, bonds, angles, dihedrals, pairs):
    # Get current atom set.
    atoms = {}
    for atid in particle_ids:
        p = system.storage.getParticle(atid)
        atom_params = topol.atom_type_params[p.type]
        topo_atom = TopoAtom()
        topo_atom.atom_id = atid
        topo_atom.atom_type = atom_params['type']
        topo_atom.chain_name = atom_params['molecule']
        topo_atom.name = 'T{}'.format(p.type)
        topo_atom.mass = atom_params['mass']
        topo_atom.charge = atom_params['charge']
        topo_atom.chain_idx = p.res_id
        topol.topol.atoms[atid] = topo_atom

    for fpl in bonds:
        for bp in fpl.getBonds():
            for b12 in bp:
                topol.topol.new_data['bonds'][b12] = []

    for ftl in angles:
        for bt in ftl.getTriples():
            for b123 in bt:
                topol.topol.angles[b123] = []

    for fql in dihedrals:
        for bq in fql.getQuadruples():
            for b1234 in bq:
                topol.topol.dihedrals[b1234] = []

    for fpr in pairs:
        for bl in fpr.getBonds():
            for b12 in bl:
                topol.topol.pairs[b12] = []

    topol.topol.write(file_name)
