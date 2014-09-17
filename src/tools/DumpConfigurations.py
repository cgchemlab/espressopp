#  Copyright (C) 2012,2013
#      Max Planck Institute for Polymer Research
#  Copyright (C) 2008,2009,2010,2011
#      Max-Planck-Institute for Polymer Research & Fraunhofer SCAI
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

import espresso
import math


def writegro(filename, system, N_mol, velocities=True, unfolded=False,
             append=False, nameseq=None, scale=1.0, t=None):
  """Writes the Gromacs .gro file.

  Args:
    filename: The file name of .gro file.
    system: The system object.
    N_mol: Number of atoms in single molecule.
    velocities: True for storing the velocities.
    unfolded: True if the position should be stored without BC.
    append: True if the file should be appended with new position.
    nameseq: The optional list of atom names to be used.
    scale: The optional value for rescaling the values.
    t: The optional value of time.
  """

  if append:
    gro_file = open(filename, 'a')
  else:
    gro_file = open(filename, 'w')

  if t is None:
    gro_file.write('MD Espresso++\n')
  else:
    gro_file.write('MD Espresso++, t= %f\n' % t)

  num_particles = int(espresso.analysis.NPart(system).compute())
  gro_file.write('%d\n' % num_particles)

  box = system.bc.boxL
  if nameseq is None:
    nameseq = ['X']
  nameseq_l = len(nameseq)

  st = '%5d%-5s%5s%5d%8.3f%8.3f%8.3f%8.4f%8.4f%8.4f\n'
  configurations = espresso.analysis.ConfigurationsExt(system)
  configurations.unfolded = unfolded
  configurations.gather()
  configuration = configurations[0]

  if velocities:
    velocities = espresso.analysis.Velocities(system)
    velocities.gather()
    velocity = velocities[0]

  for pid in configuration:
    xpos = configuration[pid][0] * scale
    ypos = configuration[pid][1] * scale
    zpos = configuration[pid][2] * scale

    if velocities:
      xvel = velocity[pid][0] * scale
      yvel = velocity[pid][1] * scale
      zvel = velocity[pid][2] * scale
    else:
      xvel = yvel = zvel = 0.0
    mol_id = int(math.ceil((pid - 1) / N_mol) + 1)
    gro_file.write(st % (
        mol_id, 'UNX', nameseq[(pid - 1) % nameseq_l], pid, xpos, ypos, zpos, xvel, yvel, zvel
        ))

  gro_file.write('%-15.10f %15.10f %15.10f\n' % (box[0], box[1], box[2]))
  gro_file.close()


def writexyz(filename, system, velocities=True, unfolded=False, append=False):
  """Writes the xyz file."""
  if append:
    input_file = open(filename, 'a')
  else:
    input_file = open(filename, 'w')
  numParticles = int(espresso.analysis.NPart(system).compute())
  box_x = system.bc.boxL[0]
  box_y = system.bc.boxL[1]
  box_z = system.bc.boxL[2]
  st = "%d\n%15.10f %15.10f %15.10f\n" % (numParticles, box_x, box_y, box_z)
  input_file.write(st)
  maxParticleID = int(espresso.analysis.MaxPID(system).compute())
  pid = 0
  while pid <= maxParticleID:
    if system.storage.particleExists(pid):
      particle = system.storage.getParticle(pid)
      if unfolded == False:
        xpos = particle.pos[0]
        ypos = particle.pos[1]
        zpos = particle.pos[2]
      else:
        unfoldedpos = system.bc.getUnfoldedPosition(particle.pos, particle.imageBox)
        xpos = unfoldedpos[0]
        ypos = unfoldedpos[1]
        zpos = unfoldedpos[2]
      xvel = particle.v[0]
      yvel = particle.v[1]
      zvel = particle.v[2]
      ptype = particle.type
      if velocities:
        st = "%d %d %15.10f %15.10f %15.10f %15.10f %15.10f %15.10f\n" % (
            pid, ptype, xpos, ypos, zpos, xvel, yvel, zvel)
      else:
        st = "%d %d %15.10f %15.10f %15.10f\n" % (pid, ptype, xpos, ypos, zpos)
      input_file.write(st)
      pid += 1
    else:
      pid += 1

  file.close()


def readxyz(filename):
  input_file = open(filename)
  line = input_file.readline()
  num_particles = int(line.split()[0])
  line = input_file.readline().split()
  if len(line) == 3:
    Lx = float(line[0])
    Ly = float(line[1])
    Lz = float(line[2])
  else:
    Lx = float(line[0])
    Ly = float(line[4])
    Lz = float(line[8])

  pid = []
  ptype = []
  xpos = []
  ypos = []
  zpos = []
  xvel = []
  yvel = []
  zvel = []
  for _ in range(num_particles):
    line = input_file.readline().split()
    if len(line) == 7 or len(line) == 4:
      line.insert(1, '0')
    pid.append(int(line[0]))
    ptype.append(int(line[1]))
    xpos.append(float(line[2]))
    ypos.append(float(line[3]))
    zpos.append(float(line[4]))
    if len(line) > 5:
      xvel.append(float(line[5]))
      yvel.append(float(line[6]))
      zvel.append(float(line[7]))
    else:
      xvel.append(0.0)
      yvel.append(0.0)
      zvel.append(0.0)

  input_file.close()
  return pid, ptype, xpos, ypos, zpos, xvel, yvel, zvel, Lx, Ly, Lz

def readxyzr(filename):
  input_file = open(filename)
  line = input_file.readline()
  num_particles = int(line.split()[0])
  line = input_file.readline()
  Lx = float(line.split()[0])
  Ly = float(line.split()[1])
  Lz = float(line.split()[2])
  pid = []
  ptype = []
  xpos = []
  ypos = []
  zpos = []
  xvel = []
  yvel = []
  zvel = []
  radius = []
  for _ in range(num_particles):
    line = input_file.readline().split()
    if len(line) == 7:
      line.insert(1, '0')
    pid.append(int(line[0]))
    ptype.append(int(line[1]))
    xpos.append(float(line[2]))
    ypos.append(float(line[3]))
    zpos.append(float(line[4]))
    if len(line) == 6:
      radius.append(float(line[5]))
    else:
      radius.append(0.0)
    if len(line) > 5 and len(line) <= 8 and len(line) != 6:
      xvel.append(float(line[5]))
      yvel.append(float(line[6]))
      zvel.append(float(line[7]))
    else:
      xvel.append(0.0)
      yvel.append(0.0)
      zvel.append(0.0)
    if len(line) == 9:
      radius.append(float(line[8]))
    else:
      if len(line) != 6:
        radius.append(0.0)
  input_file.close()
  return pid, ptype, xpos, ypos, zpos, xvel, yvel, zvel, Lx, Ly, Lz, radius


# Livia's modified writexyz to fastwritexyz with velocities
def fastwritexyz(filename, system, velocities=True, unfolded=True, append=False, scale=1.0):

  if append:
    input_file = open(filename, 'a')
  else:
    input_file = open(filename, 'w')

  configurations = espresso.analysis.ConfigurationsExt(system)
  configurations.unfolded = unfolded
  configurations.gather()
  configuration = configurations[0]

  if velocities:
    velocities = espresso.analysis.Velocities(system)
    velocities.gather()
    velocity = velocities[0]

  numParticles = int(espresso.analysis.NPart(system).compute())
  box_x = system.bc.boxL[0]*scale
  box_y = system.bc.boxL[1]*scale
  box_z = system.bc.boxL[2]*scale
  st = "%d\n%15.10f %15.10f %15.10f\n" % (numParticles, box_x, box_y, box_z)
  input_file.write(st)

  for pid in configuration:
    xpos = configuration[pid][0]*scale
    ypos = configuration[pid][1]*scale
    zpos = configuration[pid][2]*scale
    if velocities:
      xvel = velocity[pid][0]*scale
      yvel = velocity[pid][1]*scale
      zvel = velocity[pid][2]*scale
      st = "%d %15.10f %15.10f %15.10f %15.10f %15.10f %15.10f\n" % (
          pid, xpos, ypos, zpos, xvel, yvel, zvel)
    else:
      st = "%d %15.10f %15.10f %15.10f\n" % (
          pid, xpos, ypos, zpos)
    input_file.write(st)
    #pid   += 1

  input_file.close()


# Franziska's modified readxyz to fastreadxyz without velocities
def fastreadxyz(filename):
  input_file = open(filename)
  line = input_file.readline()
  num_particles = int(line.split()[0])
  line = input_file.readline().split()
  if len(line) == 3:
    Lx = float(line[0])
    Ly = float(line[1])
    Lz = float(line[2])
  else:
    Lx = float(line[0])
    Ly = float(line[4])
    Lz = float(line[8])

  pid = []
  ptype = []
  xpos = []
  ypos = []
  zpos = []
  for _ in range(num_particles):
    line = input_file.readline().split()
    if len(line) == 7 or len(line) == 4:
      line.insert(1, '0')
    pid.append(int(line[0]))
    ptype.append(int(line[1]))
    xpos.append(float(line[2]))
    ypos.append(float(line[3]))
    zpos.append(float(line[4]))
  input_file.close()
  return pid, ptype, xpos, ypos, zpos, Lx, Ly, Lz

def fastwritexyz_standard(filename, system, unfolded=False, append=False):
  """Fast write standard xyz input_file. Generally standard xyz file is
  >  number of particles
  >  comment line
  >  ptype x y z
  >  ......
  >  ......
  >  ......

  Additional information can be found here:
  Wiki:  http://en.wikipedia.org/wiki/XYZ_input_file_format
  OpenBabel: http://openbabel.org/wiki/XYZ_%28format%29

  In this case one can choose folded or unfolded coordinates.
  Currently it writes only particle ptype = 0 and pid is a line number.
  Later different ptypes should be implemented.
  """

  if append:
    input_file = open(filename, 'a')
  else:
    input_file = open(filename, 'w')

  conf = espresso.analysis.ConfigurationsExt(system)
  conf.unfolded = unfolded
  conf.gather()

  numParticles = int(espresso.analysis.NPart(system).compute())
  box_x = system.bc.boxL[0]
  box_y = system.bc.boxL[1]
  box_z = system.bc.boxL[2]
  st = "%d\n%18.12f %18.12f %18.12f\n" % (numParticles, box_x, box_y, box_z)
  input_file.write(st)

  for pid in conf[0]:
    xpos = conf[0][pid][0]
    ypos = conf[0][pid][1]
    zpos = conf[0][pid][2]

    st = "%d %15.10f %15.10f %15.10f\n"%(0, xpos, ypos, zpos)
    input_file.write(st)
    pid += 1

  input_file.close()


def xyzinput_filewrite(filename, system, append=False,
                       atomtypes=None, velocities=False, charge=False):
  """This method creates a xyz input_file with the data from a specific system:
  1. row:         number of the atoms
  2. row:         REMARK generated by ESPResSo++
  following rows: atomsymbol positionX positionY positionZ (velocityX velocityY velocityZ) (charge)
  last row:       END

  The method needs the following parameters:

  * filename
      name of the input_file where the table schould be saved in
  * system
      ESPResSo system which creates the data e.g.:
      >>>system, integrator = espresso.standard_system.LennardJones(100,(10,10,10))
  * append
      =False
        the data in the input_file will be overwritten
      =True
        the data will be appended
  * atomtypes
      the xyz input_file needs atom symbols, so it has to translate the numbers
      insert a dictionary with the right translation
  * velocities
      =False
        does not save the velocity vectors
      =True
        creates collumns for the velocity vectors and saves the data
  * charge
      =False
        does not save the charge
      =True
        creates collumns for the charges and saves the data
  """
  if atomtypes is None:
    atomtypes = {0:'Fe', 1:'O', 2:'C'}
  if append:
    input_file = open(filename, 'a')
  else:
    input_file = open(filename, 'w')
  maxParticleID = int(espresso.analysis.MaxPID(system).compute())
  pid = 0
  comment = "REMARK generated by ESPResSo++"

  st = "%d\n%s\n" % (maxParticleID, comment)
  input_file.write(st)

  while pid <= maxParticleID:
    if system.storage.particleExists(pid):
      particle = system.storage.getParticle(pid)

      xpos = particle.pos[0]
      ypos = particle.pos[1]
      zpos = particle.pos[2]
      ptype = particle.type
      vx = particle.v[0]
      vy = particle.v[1]
      vz = particle.v[2]
      q = particle.q

      if ptype in atomtypes:
        atom = atomtypes[ptype]
      else:
        atom = 'XX'

    if velocities == True:
      if charge == True:
        st = "%s %15.10f %15.10f %15.10f %15.10f %15.10f %15.10f %15.10f\n" % (
            atom, xpos, ypos, zpos, vx, vy, vz, q)
      else:
        st = "%s %15.10f %15.10f %15.10f %15.10f %15.10f %15.10f\n" % (
            atom, xpos, ypos, zpos, vx, vy, vz)
    else:
      if charge == True:
        st = "%s %15.10f %15.10f %15.10f %15.10f\n" % (atom, xpos, ypos, zpos, q)
      else:
        st = "%s %15.10f %15.10f %15.10f\n" % (atom, xpos, ypos, zpos)

    input_file.write(st)
    pid += 1

  input_file.write('END\n')
  input_file.close()

