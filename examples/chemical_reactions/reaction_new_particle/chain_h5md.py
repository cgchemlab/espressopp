# Pierre de Buyl 2014
# This file is licensed under the modified BSD license

import espressopp
from espressopp import Real3D
import pyh5md
import numpy as np
import h5py

def DumpH5MD(filename, system, integrator, author, author_email=None, edges=None, edges_time=False, n_states=None, valid_types=None):
    espressopp.Version().info()
    f = pyh5md.H5MD_File(filename, 'w', creator='espressopppp', creator_version=espressopp.Version().info(), author=author, author_email=author_email)
    if valid_types is None:
        valid_types = []
    f.valid_types = valid_types
    atoms = f.particles_group('atoms')
    maxParticleID = int(espressopp.analysis.MaxPID(system).compute())
    pos = atoms.trajectory('position', (maxParticleID+1,3), np.float64)
    species = atoms.trajectory('species', (maxParticleID+1,), np.int32)
    state = atoms.trajectory('state', (maxParticleID+1,), np.int32)
    lambda_adr = atoms.trajectory('lambda_adr', (maxParticleID+1,), np.float64)
    res_id = atoms.trajectory('res_id', (maxParticleID+1, ), np.int32)
    if edges_time:
        f.box = atoms.box(dimension=3, boundary=['periodic', 'periodic', 'periodic'], edges=edges, time=True)
    else:
        f.box = atoms.box(dimension=3, boundary=['periodic', 'periodic', 'periodic'], edges=edges)
    f.NPart = espressopp.analysis.NPart(system).compute()
    
    f.observable('particle_number', data=int(f.NPart), time=False)
    f.f['observables'].attrs['dimension']=3
    obs_dict = {}
    f.n_states = n_states
    state_tuple = (('statecount', (n_states,), np.int32),) if n_states is not None else ()
    for o in (
            ('temperature', (), np.float64),
            ('kinetic_energy', (), np.float64),
            ('pressure', (), np.float64), 
            ('pressure_tensor', (6,), np.float64),
            ('potential_energy', (), np.float64),
            ('internal_energy', (), np.float64),
            ('lennard_jones', (), np.float64),
            ('lennard_jones_dynamic', (), np.float64),
            ('harmonic', (), np.float64),
    ) + state_tuple:
        obs_dict[o[0]] = f.observable(*o)
    
    f.parameters = f.f.create_group('parameters')
   
    def dump():
        step = integrator.step
        time = integrator.step*integrator.dt
        maxParticleID = int(espressopp.analysis.MaxPID(system).compute())
        if maxParticleID>pos.value.shape[1]:
            raise ValueError('System too large for dataset')
        tmp_r = []
        tmp_species = []
        tmp_state = []
        tmp_lambda = []
        tmp_res_id = []
        for pid in range(maxParticleID+1):
            p = system.storage.getParticle(pid)
            tmp_r.append([x for x in p.pos])
            tmp_species.append(p.type)
            tmp_state.append(p.state)
            tmp_lambda.append(p.lambda_adr)
            tmp_res_id.append(p.res_id)
        pos.append(np.array(tmp_r), step, time)
        species.append(np.array(tmp_species), step, time)
        state.append(np.array(tmp_state), step, time)
        lambda_adr.append(np.array(tmp_lambda), step, time)
        res_id.append(np.array(tmp_res_id), step, time)
    f.dump = dump
    def flush():
        f.f.flush()

    f.flush = flush

    def analyse():
        step = integrator.step
        time = integrator.step*integrator.dt
        T_comp = espressopp.analysis.Temperature(system)
        for t in f.valid_types:
            T_comp.add_type(t)
        T = T_comp.compute()
        obs_dict['temperature'].append(T, step, time)
        P      = espressopp.analysis.Pressure(system).compute()
        obs_dict['pressure'].append(P, step, time)
        Pij    = espressopp.analysis.PressureTensor(system).compute()
        obs_dict['pressure_tensor'].append(Pij, step, time)
        Ek     = (3.0/2.0) * T
        obs_dict['kinetic_energy'].append(Ek, step, time)
        Harmonic = system.getInteraction(0).computeEnergy()
        obs_dict['harmonic'].append(Harmonic, step, time)
        LJ = system.getInteraction(1).computeEnergy()
        obs_dict['lennard_jones'].append(LJ, step, time)
        DynLJ = system.getInteraction(2).computeEnergy()
        obs_dict['lennard_jones_dynamic'].append(DynLJ, step, time)
        potential = LJ+DynLJ+Harmonic
        obs_dict['potential_energy'].append(potential, step, time)
        obs_dict['internal_energy'].append(Ek+potential, step, time)
        if f.n_states is not None:
            obs_dict['statecount'].append(np.bincount(np.array([system.storage.getParticle(pid).state for pid in range(maxParticleID+1)]), minlength=f.n_states), step, time)
    f.analyse = analyse
    return f

class DumpTopo(object):
    def __init__(self, f, group, name, system, integrator, fpl, time=False, chunks=None):
        self.file = f
        self.system = system
        self.integrator = integrator
        self.fpl = fpl
        f.f.require_group('topology')
        f.f.require_group('topology/'+group)
        if time:
            self.element = pyh5md.base.TimeData(f.f['topology/atoms'], name, shape=(0,2), dtype=np.int32, chunks=chunks, fillvalue=-1)
        else:
            data = np.array([b for local_bonds in fpl.getBonds() for b in local_bonds])
            self.element = pyh5md.base.FixedData(f.f['topology/atoms'], name, data=data)

    def dump(self):
        if not isinstance(self.element, pyh5md.base.TimeData):
            raise UserWarning("Trying to append data to a non suitable object.")
        bl = np.array([b for local_bonds in self.fpl.getBonds() for b in local_bonds])
        self.element.append(bl, self.integrator.step, self.integrator.step*self.integrator.dt)

