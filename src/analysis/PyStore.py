
import espressopp
from espressopp.esutil import cxxinit
from espressopp import pmi
from _espressopp import analysis_PyStore

from mpi4py import MPI
import numpy as np
import pyh5md


class PyStoreLocal(analysis_PyStore):
    def __init__(self, system, filename, group_name='all',
                 store_position=True, store_species=True,
                 store_state=False, store_velocity=False,
                 store_force=False,
                 store_charge=False,
                 store_lambda=False,
                 static_box=True,
                 is_adress=False,
                 author='xxx',
                 email='xxx',
                 chunk_size=256):
        """
        Args:
            system: The system object.
            filename: The name of hdf file name.
            group_name: The name of atom groups. (default: 'all').
            store_position: If set to True then position will be stored. (default: True)
            store_species: If set to True then species will be stored. (default: True)
            store_state: If set to True then state will be stored. (default: False)
            store_velocity: If set to True then velocity will be stored. (default: False)
            store_force: If set to True then force will be stored. (default: False)
            store_charge: If set to True then charge will be stored. (default: False)
            store_lambda: If set to True then lambda (AdResS) will be stored. (default: False)
            static_box: If set to True then box is static (like in NVT ensemble) (default: True)
            is_adress: If set to True then AdResS particles will be save instead of
                coarse-grained.
            author: The name of author of the file. (default: xxx)
            email: The e-mail to author of that file. (default: xxx)
            chunk_size: The size of data chunk. (default: 256)
        """
        if not pmi.workerIsActive():
            return
        cxxinit(self, analysis_PyStore, system, is_adress)

        self.group_name = group_name
        self.store_position = store_position
        self.store_species = store_species
        self.store_state = store_state
        self.store_velocity = store_velocity
        self.store_force = store_force
        self.store_charge = store_charge
        self.store_lambda = store_lambda
        self.static_box = static_box
        self.chunk_size = chunk_size

        self.system = system
        self.file = pyh5md.H5MD_File(filename, 'w', driver='mpio', comm=MPI.COMM_WORLD,
                                     creator='espressopp',
                                     creator_version=espressopp.VersionLocal().info(),
                                     author=author, email=email
                                     )

        self._system_data()

        part = self.file.particles_group(self.group_name)
        if self.static_box:
            self.box = part.box(dimension=3,
                                boundary=['periodic', 'periodic', 'periodic'],
                                time=False,
                                edges=np.array(
                                    [ed_i for ed_i in self.system.bc.boxL],
                                    dtype=np.float64
                                ))
        else:
            self.box = part.box(
                dimension=3,
                boundary=['periodic', 'periodic', 'periodic'],
                time=True,
                edges=np.zeros(3, dtype=np.float64))
        self.id_e = part.trajectory(
            'id', (self.chunk_size,), np.int, chunks=(1, self.chunk_size), fillvalue=-1)
        self.mass = part.trajectory(
            'mass', (self.chunk_size,), np.float64, chunks=(1, self.chunk_size), fillvalue=-1)
        if self.store_position:
            self.position = part.trajectory(
                'position', (self.chunk_size, 3), np.float64, chunks=(1, self.chunk_size, 3))
            self.image = part.trajectory(
                'image', (self.chunk_size, 3), np.float64, chunks=(1, self.chunk_size, 3))
            self.res_id = part.trajectory(
                'res_id', (self.chunk_size,), np.int, chunks=(1, self.chunk_size), fillvalue=-1)
        if self.store_species:
            self.species = part.trajectory(
                'species', (self.chunk_size,), np.int, chunks=(1, self.chunk_size), fillvalue=-1)
        if self.store_state:
            self.state = part.trajectory(
                'state', (self.chunk_size,), np.int, chunks=(1, self.chunk_size), fillvalue=-1)
        if self.store_velocity:
            self.velocity = part.trajectory(
                'velocity', (self.chunk_size, 3), np.float64, chunks=(1, self.chunk_size, 3))
        if self.store_force:
            self.force = part.trajectory(
                'force', (self.chunk_size, 3), np.float64, chunks=(1, self.chunk_size, 3))
        if self.store_charge:
            self.charge = part.trajectory(
                'charge', (self.chunk_size,), np.float64, chunks=(1, self.chunk_size), fillvalue=-1)
        if self.store_lambda:
            self.lambda_adr = part.trajectory(
                'lambda_adr', (self.chunk_size,), np.float64,
                chunks=(1, self.chunk_size), fillvalue=-1)

        self.connectivity_map = []
        self.connectivity_g = self.file.f.create_group('connectivity')

    def _system_data(self):
        """Stores specific information about simulation."""
        # Creates /system group
        sys_group = self.file.f.create_group('parameters')
        sys_group.attrs['software-id'] = 'espressopp'
        sys_group.attrs['rng-seed'] = self.system.rng.get_seed()
        sys_group.attrs['skin'] = self.system.skin
        if self.system.integrator is not None:
            sys_group.attrs['dt'] = self.system.integrator.dt

    def update(self):
        if pmi.workerIsActive():
            self.cxxclass.update(self)

    def clear_buffers(self):
        self.cxxclass.clear_buffers(self)

    def getPosition(self):
        return self.cxxclass.getPosition(self)

    def getImage(self):
        return self.cxxclass.getImage(self)

    def getResID(self):
        return self.cxxclass.getResID(self)

    def getVelocity(self):
        return self.cxxclass.getVelocity(self)

    def getForce(self):
        return self.cxxclass.getForce(self)

    def getId(self):
        return self.cxxclass.getId(self)

    def getSpecies(self):
        return self.cxxclass.getSpecies(self)

    def getState(self):
        return self.cxxclass.getState(self)

    def getCharge(self):
        return self.cxxclass.getCharge(self)

    def getMass(self):
        return self.cxxclass.getMass(self)

    def getLambdaAdr(self):
        return self.cxxclass.getLambda(self)

    def dump(self, step, time):
        if not pmi.workerIsActive():
            return
        self.update()
        NLocal = np.array(self.NLocal, 'i')
        NMaxLocal = np.array(0, 'i')
        MPI.COMM_WORLD.Allreduce(NLocal, NMaxLocal, op=MPI.MAX)
        cpu_size = ((NMaxLocal//self.chunk_size)+1)*self.chunk_size
        total_size = MPI.COMM_WORLD.size*cpu_size
        idx_0 = MPI.COMM_WORLD.rank*cpu_size
        idx_1 = idx_0+NLocal

        # Store ids.
        id_ar = np.asarray(self.getId())
        if total_size > self.id_e.value.shape[1]:
            self.id_e.value.resize(total_size, axis=1)
        self.id_e.append(id_ar, step, time, region=(idx_0, idx_1))

        if self.store_position:
            pos = np.asarray(self.getPosition())
            if total_size > self.position.value.shape[1]:
                self.position.value.resize(total_size, axis=1)
            self.position.append(pos, step, time, region=(idx_0, idx_1))
            if not self.static_box:
                self.box.edges.append(
                    np.array([edge_i for edge_i in self.system.bc.boxL], dtype=np.float64),
                    step,
                    time)
            # Store image.
            image = np.asarray(self.getImage())
            if total_size > self.image.value.shape[1]:
                self.image.value.resize(total_size, axis=1)
            self.image.append(image, step, time, region=(idx_0, idx_1))

            # Store res_id.
            res_id = np.asarray(self.getResID())
            if total_size > self.res_id.value.shape[1]:
                self.res_id.value.resize(total_size, axis=1)
            self.res_id.append(res_id, step, time, region=(idx_0, idx_1))

        # Store velocity.
        if self.store_velocity:
            vel = np.asarray(self.getVelocity())
            if total_size > self.velocity.value.shape[1]:
                self.velocity.value.resize(total_size, axis=1)
            self.velocity.append(vel, step, time, region=(idx_0, idx_1))

        if self.store_force:
            force = np.asarray(self.getForce())
            if total_size > self.force.value.shape[1]:
                self.force.value.resize(total_size, axis=1)
            self.force.append(force, step, time, region=(idx_0, idx_1))

        if self.store_charge:
            charge = np.asarray(self.getCharge())
            if total_size > self.charge.value.shape[1]:
                self.charge.value.resize(total_size, axis=1)
            self.charge.append(charge, step, time, region=(idx_0, idx_1))

        # Store mass.
        mass = np.asarray(self.getMass())
        if total_size > self.mass.value.shape[1]:
            self.mass.value.resize(total_size, axis=1)
        self.mass.append(mass, step, time, region=(idx_0, idx_1))

        # Store species.
        if self.store_species:
            species = np.asarray(self.getSpecies())
            if total_size > self.species.value.shape[1]:
                self.species.value.resize(total_size, axis=1)
            self.species.append(species, step, time, region=(idx_0, idx_1))

        # Store state.
        if self.store_state:
            state = np.asarray(self.getState())
            if total_size > self.state.value.shape[1]:
                self.state.value.resize(total_size, axis=1)
            self.state.append(state, step, time, region=(idx_0, idx_1))

        # Store lambda_adr
        if self.store_lambda:
            lambda_adr = np.asarray(self.getLambda())
            if total_size > self.lambda_adr.value.shape[1]:
                self.lambda_adr.value.resize(total_size, axis=1)
            self.lambda_adr.append(lambda_adr, step, time, region=(idx_0, idx_1))

    def add_connectivity(self, fixed_list, name, fixed_list_rank=2):
        if pmi.workerIsActive():
            g_fixed_list = pyh5md.base.TimeData(
                self.connectivity_g, name, shape=(0, fixed_list_rank),
                dtype=np.int)
            self.connectivity_map.append((name, fixed_list, g_fixed_list))

    def dump_connectivity(self):
        if pmi.workerIsActive():
            # store connectivity.
            for name, fl, g_ in self.connectivity_map:
                NLocal = np.array(fl.size(), 'i')
                NMaxLocal = np.array(0, 'i')
                MPI.COMM_WORLD.Allreduce(NLocal, NMaxLocal, op=MPI.MAX)
                cpu_size = ((NMaxLocal//self.chunk_size)+1)*self.chunk_size
                total_size = MPI.COMM_WORLD.size*cpu_size
                idx_0 = MPI.COMM_WORLD.rank*cpu_size
                idx_1 = idx_0+NLocal
                items = fl.getBonds()
                if total_size > g_.value.shape[1]:
                    g_.value.resize(total_size, axis=1)
                g_.append(items, step, time, region=(idx_0, idx_1))

    def close_file(self):
        self.file.close()

    def close(self):
        self.file.close()

    def flush(self):
        self.file.flush()

if pmi.isController:
    class PyStore(object):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls='espressopp.analysis.PyStoreLocal',
            pmicall=['update', 'getPosition', 'getId', 'getSpecies', 'getState', 'getImage',
                     'getVelocity', 'getMass', 'getCharge',
                     'close_file', 'dump', 'clear_buffers', 'flush', 'close',
                     'add_connectivity', 'dump_connectivity'],
            pmiproperty=['store_position', 'store_species', 'store_state', 'store_velocity',
                         'store_charge']
        )
