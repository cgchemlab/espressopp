#! /usr/bin/env python
#
# Copyright (c) 2016 Jakub Krajniak <jkrajniak@gmail.com>
#
# Distributed under terms of the GNU GPLv3 license.
#

import os
import files_io


class FileBuffer():
    def __init__(self):
        self.linecount = 0
        self.lines = []
        self.pos = 0

    def appendline(self, line):
        self.lines.append(line)

    def readline(self):
        try:
            line = self.lines[self.pos]
        except:
            return ''
        self.pos += 1
        return line

    def readlastline(self):
        try:
            line = self.lines[self.pos-1]
        except:
            return ''
        return line

    def seek(self, p):
        self.pos = p

    def tell(self):
        return self.pos


def FillFileBuffer(fname, filebuffer, cwd=None, defines=None):
    if cwd is None:
        cwd = '.'
    if defines is None:
        defines = {}
    f = open(os.path.join(cwd, fname), 'r')
    for line in f:
        if line.startswith(';'):
            continue
        if "include" in line:
            name = line.split()[1].strip('\"')
            cwd_name = os.path.dirname(name)
            if cwd_name != '':
                cwd = cwd_name
            FillFileBuffer(name, filebuffer, cwd, defines)
        elif 'define' in line:
            t = line.strip().split()
            if len(t) > 2:
                defines[t[1]] = ' '.join(t[2:])
        else:
            l = line.rstrip('\n')
            if l:
                filebuffer.appendline(l)

    f.close()
    return


def PostProcessFileBuffer(filebuffer, defines):
    """Replace all defines with the value from the dictionary."""
    ret_fb = FileBuffer()
    define_keys = set(defines)
    for line in filebuffer.lines:
        line = line.strip()
        if line:
            if not (line.startswith(';') or line.startswith('#define')
                    or line.startswith('#include') or line.startswith('#ifdef')
                    or line.startswith('#ifndef')):
                def_key = set.intersection(set(map(str.strip, line.split())), define_keys)
                if def_key:
                    def_key = def_key.pop()
                    ret_fb.appendline(
                        line.replace(def_key, defines[def_key]))
                else:
                    ret_fb.appendline(line)
            else:
                ret_fb.appendline(line)
    return ret_fb


def convertc6c12(c6, c12, cr):
    if cr == 1:
        if c12 == 0.0:
            return 1.0, 0.0
        sig = pow(c12/c6, 1.0/6.)
        if sig > 0.0:
            eps = 0.25*c6*pow(sig, -6.0)
        else:
            eps = 0.0
        return sig, eps
    else:
        return c6, c12


class GromacsTopology:
    def __init__(self, input_topol):
        self.input_file = input_topol
        self.content = None
        self.data = {}

        self.atomsym_atomtype = {}
        self.atomtype_atomsym = {}

        self.bonds = {}
        self.angles = {}
        self.dihedrals = {}

        self.bondparams = {}
        self.angleparams = {}
        self.dihedralparams = {}

    def read(self):
        fb = FileBuffer()
        defines = {}
        FillFileBuffer(self.input_file, fb, defines=defines)
        f = PostProcessFileBuffer(fb, defines)
        self.gt = files_io.GROMACSTopologyFile(self.input_file)
        self.topol = self.gt
        self.gt.content = f.lines
        self.gt.read()
        if len(self.gt.molecules) > 1:
            raise RuntimeError('Multiple molecules are not supported')

        self._prepare_data()

    def _prepare_data(self):
        # Generate atom types from atom symbols.
        self.atomsym_atomtype = {
            k: idx for idx, k in enumerate(sorted(self.gt.atomtypes))
        }
        self.atomtype_atomsym = {v: k for k, v in self.atomsym_atomtype.items()}
        self.atomparams = {}
        self.atom_id_params = {}
        self.atom_type_params = {}
        combinationrule = self.topol.defaults['combinationrule']
        for at_id in sorted(self.gt.atoms):
            at_data = self.gt.atoms[at_id]
            at_type = self.gt.atomtypes[at_data.atom_type]
            at_key = '{}-{}'.format(at_data.chain_name, at_data.name)
            self.atomparams[at_key] = {
                'molecule': at_data.chain_name,
                'type': at_data.atom_type,
                'sig': at_type['sigma'],
                'eps': at_type['epsilon'],
                'type_id': self.atomsym_atomtype[at_data.atom_type],
                'state': at_type.get('state', 0)
            }
            if at_data.charge:
                self.atomparams[at_key]['charge'] = at_data.charge
            else:
                self.atomparams[at_key]['charge'] = at_type['charge']
            if at_data.mass:
                self.atomparams[at_key]['mass'] = at_data.mass
            else:
                self.atomparams[at_key]['mass'] = at_type['mass']
            sig, eps = convertc6c12(
                at_type['sigma'], at_type['epsilon'], combinationrule)
            self.atomparams[at_key]['sig'] = sig
            self.atomparams[at_key]['eps'] = eps
            self.atom_id_params[at_id] = self.atomparams[at_key]
            self.atom_type_params[self.atomparams[at_key]['type_id']] = self.atomparams[at_key]

        # Update non_bonded params
        for k, v in self.topol.nonbond_params.items():
            if v['func'] == 1 and self.topol.defaults['combinationrule'] == 1:
                c6 = v['params'][0]
                c12 = v['params'][1]
                sig, eps = convertc6c12(c6, c12, combinationrule)
                v['params'][0] = sig
                v['params'][1] = eps

        self._prepare_bondedparams()
        self._prepare_bondedlists()
        self._prepare_exclusionlists()

    def _prepare_bondedlists(self):
        """Replicate bonded lists."""
        n_atoms = len(self.gt.atoms)
        n_mols = self.gt.molecules.values()[0]

        self.bonds = self._replicate_lists(
            n_mols, n_atoms, sorted(self.gt.bonds))
        self.angles = self._replicate_lists(
            n_mols, n_atoms, sorted(self.gt.angles))
        self.dihedrals = self._replicate_lists(
            n_mols, n_atoms, sorted(self.gt.dihedrals))
        self.pairs = self._replicate_lists(
            n_mols, n_atoms, sorted(self.gt.pairs))

    def _prepare_exclusionlists(self):
        self.exclusions = self.bonds[:]
        self.exclusions.extend(
            [(x[0], x[2]) for x in self.angles])
        self.exclusions.extend(
            [(x[0], x[3]) for x in self.dihedrals])

    def _prepare_bondedparams(self):
        """Prepares bonded params to use with FixedListTypes interaction."""
        bcount = 0
        for i in self.gt.bondtypes:
            for j in self.gt.bondtypes[i]:
                bcount += 1
                params = self.gt.bondtypes[i][j]
                t1 = self.atomsym_atomtype[i]
                t2 = self.atomsym_atomtype[j]
                self.bondparams[(t1, t2)] = params
        assert bcount == len(self.bondparams)
        acount = 0
        for i in self.gt.angletypes:
            for j in self.gt.angletypes[i]:
                for k in self.gt.angletypes[i][j]:
                    acount += 1
                    params = self.gt.angletypes[i][j][k]
                    t1 = self.atomsym_atomtype[i]
                    t2 = self.atomsym_atomtype[j]
                    t3 = self.atomsym_atomtype[k]
                    self.angleparams[(t1, t2, t3)] = params
        assert acount == len(self.angleparams)
        dcount = 0
        for i in self.gt.dihedraltypes:
            for j in self.gt.dihedraltypes[i]:
                for k in self.gt.dihedraltypes[i][j]:
                    for l in self.gt.dihedraltypes[i][j][k]:
                        dcount += 1
                        params = self.gt.dihedraltypes[i][j][k][l]
                        t1 = self.atomsym_atomtype[i]
                        t2 = self.atomsym_atomtype[j]
                        t3 = self.atomsym_atomtype[k]
                        t4 = self.atomsym_atomtype[l]
                        self.dihedralparams[(t1, t2, t3, t4)] = params
        assert dcount == len(self.dihedralparams)

    def _replicate_lists(self, n_mols, n_atoms, input_list, shift=0):
        return [
            map(lambda x: shift+x+(mol*n_atoms), l)
            for mol in range(n_mols) for l in input_list
        ]


if __name__ == '__main__':
    gt = GromacsTopology('topol.top')
    gt.read()
    gt.topol.molecules[gt.topol.molecules.keys()[0]] = 1
    gt.topol.bonds = {}
    gt.topol.new_data['bonds'] = {(1,2): [1,2,3], (2,3): [2,3,4]}
    gt.topol.write('abc.top')

    # conf = files_io.GROFile('conf.gro')
    # conf.read()
    # import tools_sim
    # import espressopp
    # import mpi4py.MPI as MPI
    # s = espressopp.System()
    # s.rng = espressopp.esutil.RNG()
    # s.bc = espressopp.bc.OrthorhombicBC(s.rng, (10, 10, 10))
    # nodeGrid = espressopp.tools.decomp.nodeGrid(MPI.COMM_WORLD.size)
    # cellGrid = espressopp.tools.decomp.cellGrid((10, 10, 10), nodeGrid, 1.4, 0.2)
    # s.storage = espressopp.storage.DomainDecomposition(s, nodeGrid, cellGrid)
    # tools_sim.setNonbondedInteractions(s, gt, None, 1.4)
    # tools_sim.setBondInteractions(s, gt)
    # tools_sim.setAngleInteractions(s, gt)
    # tools_sim.setDihedralInteractions(s, gt)
