# Pierre de Buyl 2014
# This file is licensed under the modified BSD license

import espressopp

import numpy as np  # NOQA
import pyh5md


def DumpH5MD(filename, system, integrator, author, particle_ids,
             author_email=None, edges=None, edges_time=False, n_states=None):
    espressopp.Version().info()
    N = len(particle_ids)
    f = pyh5md.H5MD_File(
        filename,
        'w',
        creator='espressopppp',
        creator_version=espressopp.Version().info(),
        author=author,
        author_email=author_email)
    atoms = f.particles_group('atoms')
    pos = atoms.trajectory('position', (N, 3), np.float64)
    atoms.trajectory(
        'species',
        (N,),
        np.int32,
        data=np.array([system.storage.getParticle(pid).type for pid in particle_ids]),
        time=False)
    if edges_time:
        f.box = atoms.box(
            dimension=3,
            boundary=['periodic', 'periodic', 'periodic'],
            edges=edges,
            time=True)
    else:
        f.box = atoms.box(
            dimension=3,
            boundary=['periodic', 'periodic', 'periodic'],
            edges=edges)
    f.NPart = espressopp.analysis.NPart(system).compute()
    f.observable('particle_number', data=int(f.NPart), time=False)
    f.f['observables'].attrs['dimension'] = 3
    obs_dict = {}
    f.n_states = n_states
    state_tuple = ()
    if n_states is not None:
        state_tuple = (('statecount', (n_states,), np.int32),)
    for o in (
            ('temperature', (), np.float64),
            ('kinetic_energy', (), np.float64),
            ('pressure', (), np.float64),
            ('pressure_tensor', (6,), np.float64),
            ('potential_energy', (), np.float64),
            ('total_energy', (), np.float64),
    ) + state_tuple:
        obs_dict[o[0]] = f.observable(*o)

    for key in system.getInteractionLabels():
        obs_dict[key] = f.observable(key, (), np.float64)

    def dump():
        step = integrator.step
        time = integrator.step*integrator.dt
        maxParticleID = int(espressopp.analysis.MaxPID(system).compute())
        if maxParticleID > pos.value.shape[1]:
            raise ValueError('System too large for dataset')
        pos.append(
            [[x for x in system.storage.getParticle(pid).pos] for pid in particle_ids],
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
            obs_dict[interaction_label].append(energy_value, step, time)
            potential_energy += energy_value

        obs_dict['potential_energy'].append(potential_energy, step, time)
        obs_dict['total_energy'].append(Ek + potential_energy, step, time)

    f.analyse = analyse
    return f
