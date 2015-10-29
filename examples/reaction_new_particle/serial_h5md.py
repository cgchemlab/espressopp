# Pierre de Buyl 2014
# This file is licensed under the modified BSD license

import espressopp
import pyh5md
import numpy as np


def DumpH5MD(filename, system, integrator, author, author_email=None,
             edges=None, edges_time=False, n_states=None):
    espressopp.Version().info()
    f = pyh5md.H5MD_File(
        filename,
        'w',
        creator='espressopppp',
        creator_version=espressopp.Version().info(),
        author=author,
        author_email=author_email)
    atoms = f.particles_group('atoms')
    maxParticleID = int(espressopp.analysis.MaxPID(system).compute())
    pos = atoms.trajectory('position', (maxParticleID+1, 3), np.float64)
    species = atoms.trajectory('species', (maxParticleID+1,), np.int32)
    state = atoms.trajectory('state', (maxParticleID+1,), np.int32)
    if edges_time:
        f.box = atoms.box(dimension=3, boundary=['periodic', 'periodic', 'periodic'],
                          edges=edges, time=True)
    else:
        f.box = atoms.box(dimension=3, boundary=['periodic', 'periodic', 'periodic'], edges=edges)
    f.NPart = espressopp.analysis.NPart(system).compute()
    f.observable('particle_number', data=int(f.NPart), time=False)
    f.f['observables'].attrs['dimension'] = 3
    obs_dict = {}
    f.n_states = n_states
    state_tuple = (('statecount', (n_states,), np.int32),) if n_states is not None else ()
    for o in (
            ('temperature', (), np.float64), ('kinetic_energy', (), np.float64),
            ('pressure', (), np.float64), ('pressure_tensor', (6,), np.float64),
            ('potential_energy', (), np.float64), ('internal_energy', (), np.float64),
            ('total_energy', (), np.float64)
    ) + state_tuple:
        obs_dict[o[0]] = f.observable(*o)

    # Attache interactions
    for key in system.getInteractionLabels():
        obs_dict[key] = f.observable(key, (), np.float64)

    def dump():
        step = integrator.step
        time = integrator.step*integrator.dt
        maxParticleID = int(espressopp.analysis.MaxPID(system).compute())
        if maxParticleID > pos.value.shape[1]:
            raise ValueError('System too large for dataset')
        r = np.array(
            [[x for x in system.storage.getParticle(pid).pos] for pid in range(maxParticleID+1)]
            )
        pos.append(r, step, time)
        del r
        species.append(
            np.array([system.storage.getParticle(pid).type for pid in range(maxParticleID+1)]),
            step, time)
        state.append(
            np.array([system.storage.getParticle(pid).state for pid in range(maxParticleID+1)]),
            step, time)
    f.dump = dump

    def analyse():
        step = integrator.step
        time = integrator.step*integrator.dt
        T = espressopp.analysis.Temperature(system).compute() / system.kb
        obs_dict['temperature'].append(T, step, time)
        P = espressopp.analysis.Pressure(system).compute()
        obs_dict['pressure'].append(P, step, time)
        Pij = espressopp.analysis.PressureTensor(system).compute()
        obs_dict['pressure_tensor'].append(Pij, step, time)

        Ek = (3.0/2.0) * T
        obs_dict['kinetic_energy'].append(Ek, step, time)

        # Gets the potential energy
        potential_energy = 0.0
        for interaction_label, interaction_id in system.getInteractionLabels().iteritems():
            energy_value = system.getInteraction(interaction_id).computeEnergy()
            if interaction_label not in obs_dict:
                obs_dict[interaction_label] = f.observable(interaction_label, (), np.float64)
            obs_dict[interaction_label].append(energy_value, step, time)
            potential_energy += energy_value

        obs_dict['potential_energy'].append(potential_energy, step, time)
        obs_dict['total_energy'].append(Ek + potential_energy, step, time)

    f.analyse = analyse
    return f


def save_bonds(h5mdfile, bond_lists):
    h5mdfile.f.require_group('topology')
    pyh5md.base.FixedData(h5mdfile.f['topology/'], 'bonds', data=bond_lists)


class DumpTopo(object):
    def __init__(self, f, system, integrator, fpl,
                 particle_ids, time=False, chunks=None, append=False):
        self.file = f
        self.system = system
        self.integrator = integrator
        self.fpl = fpl
        self.particle_ids = particle_ids
        group = 'topology/atoms'
        f.f.require_group(group)
        self.element = pyh5md.base.TimeData(
            f.f[group],
            'crosslinks',
            shape=(0, 2),
            dtype=np.int32,
            chunks=chunks,
            compression='gzip',
            fillvalue=-1)
        self.species = pyh5md.base.TimeData(
            f.f[group],
            'species',
            shape=(len(particle_ids),),
            dtype=np.int32,
            chunks=chunks,
            compression='gzip',
            fillvalue=-1
        )
        self.species.append(
            np.array([system.storage.getParticle(x).type for x in particle_ids]),
            self.integrator.step,
            self.integrator.step*self.integrator.dt
            )
        self.res_ids = pyh5md.base.TimeData(
            f.f[group],
            'res_ids',
            shape=(len(particle_ids), ),
            dtype=np.int32,
            chunks=chunks,
            compression='gzip',
            fillvalue=-1
        )
        self.res_ids.append(
            np.array([system.storage.getParticle(x).res_id for x in particle_ids]),
            self.integrator.step,
            self.integrator.step*self.integrator.dt)

    def dump(self):
        if not isinstance(self.element, pyh5md.base.TimeData):
            raise UserWarning("Trying to append data to a non suitable object.")
        bl = np.array([b for local_bonds in self.fpl.getBonds() for b in local_bonds])
        self.element.append(bl, self.integrator.step, self.integrator.step*self.integrator.dt)
        self.species.append(
            np.array([self.system.storage.getParticle(x).type for x in self.particle_ids]),
            self.integrator.step,
            self.integrator.step*self.integrator.dt)
        self.res_ids.append(
            np.array([self.system.storage.getParticle(x).res_id for x in self.particle_ids]),
            self.integrator.step,
            self.integrator.step*self.integrator.dt)
