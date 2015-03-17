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


# -*- coding: utf-8 -*-
import math
import espresso
from topology_helper import *
from operator import itemgetter # for sorting a dict


"""This Python module allows one to use GROMACS data files as the
   input to an ESPResSo++ simulation, set interactions for given
   particle types and convert GROMACS potential tables into
   ESPResSo++ tables.
   It containts functions: read(), setInteractions(), convertTable()
   """
   
def read(gro_file, top_file="", doRegularExcl=True):
    """ Read GROMACS data files.

    Keyword arguments:
    gro_file -- contains coordinates of all particles, the number of particles, velocities and box size.
    top_file -- contains topology information. Included topology files (.itp) are also read
    doRegularExcl -- if True, exclusions are generated automatically based on the nregxcl parameter (see gromacs manual)
    """

    # read gro file
    if gro_file != "":
        f = open(gro_file)
        f.readline() # skip comment line
        total_num_particles = int(f.readline())
        
        # store coordinates and velocities
        x, y, z = [], [], []
        vx, vy, vz = [], [], []
        for i in range(total_num_particles):
            s = f.readline()[20:69]
            # coordinates
            x.append(float(s[0:8]))
            y.append(float(s[8:16]))
            z.append(float(s[16:24]))
            
            if len(s.split()) > 3:
                # velocities
                vx.append(float(s[24:32]))
                vy.append(float(s[32:40]))
                vz.append(float(s[40:49]))
        
        # store box size
        Lx, Ly, Lz = map(float, f.readline().split()) # read last line, convert to float
        f.close()


    # read top and itp files
    masses, charges = [], [] # mases and charges of the whole configuration
    types=[] # tuple: atomindex(int) to atomtypeid(int)
    bonds={} # dict: key bondtypeid value: tuple of bond pairs
    angles={} # dict: key angletype value: tuple of triples
    dihedrals = {} #same...
    exclusions=[] #list of atom pairs no considered in non-bonded interactions
    
    defaults={} # gromacs default values
    atomtypeparams={} # a dict: key atomtypeid , value : class storing actual parameters of each type e.g. c6, c12, etc..      
    bondtypeparams={} # same for bonds
    angletypeparams={} # same for angles
    dihedraltypeparams={} # same for dihedrals
    
    if top_file != "":
        #f = open(top_file)
        # FileBuffer: a class which behaves like a file, but all lines are in memory
        # we use this for emulating a 'preprocessor' which handles the #include
        # statements in the .top and .itp files
        f=FileBuffer() 
        
        FillFileBuffer(top_file, f)
                
        print "Reading top file: "+top_file
        line = ''
        itp_files = []
        a = 0
        bondtypecount, angletypecount, dihedraltypecount=0,0,0
        readdefaults, readattypes, readbdtypes, readantypes, readdhtypes  = False, False, False, False, False
        defaults={} # gromacs default values
        atomtypes={} # a dict: key atomtypename(str) value: atomtypeid(int)
        bondtypes={} # a dict: key atomindex(int),atomindex(int)  value: bondtypeid(int)
        angletypes={} # a dict: key atomindex(int), atomindex(int),atomindex(int) value: angletypeid(int)
        dihedraltypes={} # a dict: key atomtindex(int), atomindex(int), atomindex(int),atomindex(int) value: dihedraltypeid(int)
        
        # it was moved out of "if" statement
        #atomtypeparams={} # a dict: key atomtypeid , value : class storing actual parameters of each type e.g. c6, c12, etc..      
        #bondtypeparams={} # same for bonds
        #angletypeparams={} # same for angles
        #dihedraltypeparams={} # same for dihedrals
        
        #atomparams={} # key: atomindex(int) value: per atom parameters e.g. q, mass
        molecules=[]
        #molecules = {} # key: moleculeid value: name (string)
        readmolecules = False
        
        for line in f.lines:
            if line[0] == ";":  # skip comment line
                continue
                
            if 'defaults' in line: # store some gromacs default values
                readdefaults =True
                continue
                
            if readdefaults:
                if line.strip() == "" or '[' in line: # end of defaults section
                    readdefaults=False
                    print "Defaults: ", defaults
                else:
                    fields=line.split()
                    if len(fields)==5: 
                        defaults={"nbtype":fields[0], "combinationrule":fields[1],
                        "genpairs":fields[2], "fudgeLJ":fields[3], "fudgeQQ":fields[4]}
                    else: 
                        defaults={"nbtype":fields[0], "combinationrule":fields[1]} 
            
            if 'atomtypes' in line: # map atom types (espresso++ uses ints)
                readattypes = True
                print "Reading atomtypes (GROMACS: ESPResSo++): "
                continue
            
            if readattypes:
                if line.strip() == "" or '[' in line: # end of atomtypes section
                    readattypes = False
                    # prints gromacs type and esp++ type
                    for t in sorted(atomtypes.items(), key=itemgetter(1)):
                        print " %s: %d"%(t[0],t[1])
                else:
                    fields=line.split()
                    attypename = fields[0]
                
                #make a map containing the properties
                # sig, eps may be c6 and c12: this is specified in the defaults
                # and converted later
                    if len(fields)==7:
                        tmpprop={"atnum":int(fields[1]), "mass": float(fields[2]),
                        "charge":float(fields[3]), "particletype":fields[4],
                        "sig":float(fields[5]), "eps":float(fields[6])}
                    else:
                        tmpprop={"mass":float(fields[1]),
                        "charge":float(fields[2]), "particletype":fields[3],
                        "sig":float(fields[4]), "eps":float(fields[5])}
                
                    
                if attypename not in atomtypes:
                    atomtypes.update({attypename:a}) # atomtypes is used when reading the "atoms" section
                    atomtypeparams.update({a:tmpprop})
                    a += 1
                    
            if 'bondtypes' in line:
                readbdtypes = True
                continue
            
            if readbdtypes:
                if line.strip() == "" or '[' in line: # end of bondtypes section
                    readbdtypes = False
                else:
                    tmp = line.split() 
                     # i: i-atomname  i, j: j-atomname
                    i, j = atomtypes[tmp[0]], atomtypes[tmp[1]]
                    p=ParseBondTypeParam(line)                
                    #check if this type has been defined before
                    bdtypeid=FindType(p, bondtypeparams)
                    if bdtypeid==None:
                        bdtypeid=len(bondtypeparams)
                        bondtypeparams.update({bdtypeid:p})
                    
                    if i in bondtypes:
                        bondtypes[i].update({j:bdtypeid})
                    else:
                        bondtypes.update({i:{j:bdtypeid}})
            
            if 'angletypes' in line:
                readantypes = True
                continue
            
            if readantypes:
                if line.strip() == "" or '[' in line: # end of angletypes section
                    readantypes = False
                else:
                    tmp = line.split()
                    i, j, k= atomtypes[tmp[0]], atomtypes[tmp[1]], atomtypes[tmp[2]]
                    p=ParseAngleTypeParam(line)
                    
                    atypeid=FindType(p, angletypeparams)
                    if atypeid==None:
                        atypeid=len(angletypeparams)
                        angletypeparams.update({atypeid:p})
                    
                    if i in angletypes:
                        if j in angletypes[i]:
                            angletypes[i][j].update({k:atypeid})
                        else:
                            angletypes[i].update({j:{k:atypeid}})
                    else:
                        angletypes.update({i:{j:{k:atypeid}}})
                    #print "FOUND angletype: ", angletypecount, " : ", p.parameters 
                    #angletypeparams.update({angletypecount:p})
                    #angletypecount+=1
            
            if 'dihedraltypes' in line:
                readdhtypes = True
                continue
            
            if readdhtypes:
                if line.strip() == "" or '[' in line: # end of angletypes section
                    readdhtypes = False
                else:
                    tmp = line.split()
                    i, j, k, l = atomtypes[tmp[0]], atomtypes[tmp[1]], atomtypes[tmp[2]], atomtypes[tmp[3]]
                    p=ParseDihedralTypeParam(line)
                    
                    dtypeid=FindType(p, dihedraltypeparams)
                    if dtypeid==None:
                        dtypeid=len(dihedraltypeparams)
                        dihedraltypeparams.update({dtypeid:p})
                    if i in dihedraltypes:
                        if j in dihedraltypes[i]:
                            if k in dihedraltypes[i][j]:
                                dihedraltypes[i][j][k].update({l:dtypeid})
                            else:
                                dihedraltypes[i][j].update({k:{l:dtypeid}})
                        else:
                            dihedraltypes[i].update({j:{k:{l:dtypeid}}})
                    else:
                        dihedraltypes.update({i:{j:{k:{l:dtypeid}}}})
            
            if 'include' in line: # add included topology files
                itp_files.append((line.split()[1]).strip('\"')) # strip " and add filename
            
            if 'molecules' in line: # store number of chains
                readmolecules = True
                print "Reading number of molecules: "
                continue
            
            if readmolecules:
                if line.strip() == "" or '[' in line: # end of molecules section
                    readmolecules = False
                else:
                    print " "+line.strip('\n')
                    mol, nrmol = line.split()
                    #we have to check if the same molecules comes multiple times in the molecules section
                    if len(molecules) == 0:
                      molecules.append({'name':mol, 'count':int(nrmol)})
                    elif molecules[-1]['name'] == mol: #check if mol was added earlier already
                        molecules[-1]['count'] = molecules[-1]['count'] + int(nrmol) #update count
                    else: molecules.append({'name':mol, 'count':int(nrmol)}) #if mol newly added
              
        
        molstartindex=0 #this is the index of the first atom in the molecule being parsed
        
        
        f.seek(0) # Now we search for bonds, angles definitions and start from the beginning of the file buffer

        for mol in molecules:
            print "Preparing %d %s molecules... " %(mol['count'], mol['name']) 
            print "-----------------------------"

            # find and store number of molecules
            
            num_molecule_copies=mol['count']
            # this does not what the name suggests....
            nrexcl = storeMolecules(f, molecules, mol)
            # find and store atom types
            types,masses,charges, num_atoms_molecule = \
            storeAtoms(f, types, atomtypes, atomtypeparams, masses, charges, num_molecule_copies)

            # find and store bonds
            bonds = storeBonds(f, types, bondtypes, bondtypeparams, bonds,
                               num_atoms_molecule, num_molecule_copies, molstartindex, exclusions, nrexcl, doRegularExcl)

            # find and store angles
            angles = storeAngles(f, types, angletypes, angletypeparams, angles,
                                 num_atoms_molecule, num_molecule_copies, molstartindex)

            # find and store dihedrals
            dihedrals = storeDihedrals(f, types, dihedraltypes, dihedraltypeparams, dihedrals,
                                       num_atoms_molecule, num_molecule_copies, molstartindex)

            molstartindex+=num_molecule_copies*num_atoms_molecule
            
            
    
    params = []
    
    unpackvars=[]
    
    # The data is packed into a touple, unpackvars contains a string which
    # tells the user which kind of data was read.
    
    if len(defaults) != 0:
        print "Found default values"
        unpackvars.append("defaults")
        params.append(defaults)
    if len(types) != 0:
        print "Found ", len(types), "types"
        unpackvars.append("types")
        params.append(types)
    if len(masses) != 0:
        print "Found ", len(masses), "masses"
        unpackvars.append("masses")
        params.append(masses)
    if len(charges) != 0:
        print "Found ", len(charges), "charges"
        unpackvars.append("charges")
        params.append(charges)
    if len(atomtypeparams) !=0:
        print "Found ", len(atomtypeparams), " atomtypeparameters" 
        unpackvars.append("atomtypeparameters")
        params.append(atomtypeparams)
    if len(bonds) != 0:
        print "Found ", len(bonds), " bond types"
        unpackvars.append("bondtypes")
        params.append(bonds)
    if len(bondtypeparams) !=0:
        print "Found ", len(bondtypeparams), " bondtypeparams" 
        unpackvars.append("bondtypeparams")
        params.append(bondtypeparams)        
    if len(angles) != 0:
        print "Found ", len(angles), " angle types"
        unpackvars.append("angletypes")
        params.append(angles)
    if len(angletypeparams) != 0:
        print "Found ", len(angletypeparams), " angle type parameters"
        unpackvars.append("angletypeparams")
        params.append(angletypeparams)      
    if len(dihedrals) != 0:
        unpackvars.append("dihedraltypes")
        print "Found ", len(dihedrals), " dihedral types"
        params.append(dihedrals)
    if len(dihedraltypeparams) != 0:
        print "Found ", len(dihedraltypeparams), " dihedral type parameters"
        unpackvars.append("dihedraltypeparams")
        params.append(dihedraltypeparams)       
    if len(exclusions) != 0:
        print "Found ", len(exclusions), "bond exclusions"
        unpackvars.append("exclusions")
        params.append(exclusions)  
        
    unpackvars.append("x, y, z")
    params.extend([x, y, z])
    print "Found Box:", [Lx, Ly, Lz]
    if len(vx) != 0:
        params.extend([vx, vy, vz])
        print "Found ", len(vx), " velocities"
        unpackvars.append("vx, vy, vz")
        
    params.extend([Lx, Ly, Lz])
    unpackvars.append("Lx, Ly, Lz")
    
    print "USAGE: unpack as"
    s=""
    for i in range(len(unpackvars)):
        s+=str(unpackvars[i])
        if (i< len(unpackvars)-1): s+=", "
    print s, "=gromacs.read( ... )"
    return tuple(params)


def storeMolecules(f, molecules, mol=""):
    nrexcl=0
    line = ''
    line=f.readlastline()
    while not 'moleculetype' in line:
        line = f.readline()
        if not line: break # break out of while if EOF
    line = f.readline()
    while(len(line) > 1 and not '[' in line):
        if line[0] == ";":   # skip comment lines
            #print "skipping line: "+line.strip("\n")
            line = f.readline()
            continue
        fields=line.split()
        #mol = fields[0]
        nrexcl=int(fields[1])
        line = f.readline()
    return nrexcl


def storeAtoms(f, types, atomtypes, atomtypeparams, masses, charges, num_molecule_copies):
    line = ''
    types_tmp = []
    charge_tmp =[]
    mass_tmp=[]
    
    line=f.readlastline()
    while not 'atoms' in line:
        line = f.readline()
        if not line: break # break out of while if EOF
    line = f.readline()
    while(len(line) > 1 and not '[' in line):
        if line[0] == ";":   # skip comment lines
            line = f.readline()
            continue
        fields=line.split()
        attypeid=atomtypes[fields[1]] # map str type to int type
        types_tmp.append(attypeid) 
        if len(fields) > 6:
            # this atom has a charge different from its atomtype
            charge_tmp.append(float(fields[6]))
        else:
            #look up default values for this atom type
            charge_tmp.append(atomtypeparams[attypeid]['charge'])
        if len(fields) > 7:
            # also has a special mass
            mass_tmp.append(float(fields[7]))
        else:
            mass_tmp.append(atomtypeparams[attypeid]['mass'])
            
        line = f.readline()
    
    # extend copies of this molecule
    num_atoms_molecule = len(types_tmp)
    for i in range(num_molecule_copies):
        types.extend(types_tmp)
        charges.extend(charge_tmp)
        masses.extend(mass_tmp)
    
    return types, masses, charges, num_atoms_molecule

def storeBonds(f, types, bondtypes, bondtypeparams, bonds, num_atoms_molecule,\
    num_molecule_copies, molstartindex, exclusions, nregxcl, doRegularExcl=True):
    line = ''
    bonds_tmp = []
    top = False
    pos = f.tell()
    line=f.readlastline()
    local_exclusions=[] # excluded pairs of atoms within this mol (local ids)
      
    while not 'bonds' in line:
        line = f.readline()
        if 'moleculetype' in line or not line:
            f.seek(pos)
            return bonds

    line = f.readline()
    while(len(line) > 1 and not '[' in line):
        if line[0] == ";":   # skip comment lines
            line = f.readline()
            continue
        tmp = line.split()
        lookup=(len(tmp)<=3) # if the bond has < 3 arguments, it is defined in the bondtypes section and we have to look it up
        pid1, pid2 = map(int, tmp[0:2])
        if lookup:
            # based on atom names: potential has to be defined in bondtypes already
            # this is for tabulated bond potentials specified based on type
            t1, t2 = types[pid1-1], types[pid2-1]
            if t1 > t2: # interactions in the other way
                t1, t2 = t2, t1
            bdtypeid = bondtypes[t1][t2] #bondtypes[t1][t2]
        else:
            # this one is specific for this pair of atoms: check if we need to make a new type
            temptype=ParseBondTypeParam(line)
            bdtypeid=FindType(temptype, bondtypeparams)
            if bdtypeid==None:
                bdtypeid=len(bondtypeparams)
                bondtypeparams.update({bdtypeid:temptype})
            
        bonds_tmp.append((pid1, pid2, bdtypeid)) # store bondtypes for this molecule
        if bondtypeparams[bdtypeid].automaticExclusion():
             # this bond type generates an exclusion as defined by the
             # function type (see gromacs manual)
            local_exclusions.append((pid1, pid2))
        line = f.readline()
    

    if doRegularExcl:
        # generate exclusions for atoms up to a number of nregxcl bonds away
        # see gromacs manual, section 5.4
        exclusions_bonds=[]
        for b in bonds_tmp:
            pid1, pid2, bdtypeid = b[0:3]
            exclusions_bonds.append((pid1, pid2))   
        print "Generating Regular exclusions nregxcl=", nregxcl
        local_exclusions=GenerateRegularExclusions(exclusions_bonds, nregxcl,local_exclusions)
    # extend bonds to copies of this molecule
    bonds_per_mol = len(bonds_tmp)
    for i in range(num_molecule_copies):
        for j in range(bonds_per_mol):
            pid1, pid2, bdtypeid = bonds_tmp[j][0:3]
            ia=molstartindex+pid1 + (i * num_atoms_molecule) # index of copy atom i
            ib=molstartindex+pid2 + (i * num_atoms_molecule) # index of copy atom j
            
            if bdtypeid in bonds:
                bonds[bdtypeid].append((ia, ib))
            else:
                bonds.update({bdtypeid:[(ia, ib)]})
                
                
    # now, extend also the regular exclusions
    for i in range(num_molecule_copies):
        for exclpair in local_exclusions:
            pid1, pid2 = exclpair[0:2]
            ia=molstartindex+pid1 + (i * num_atoms_molecule) # index of copy atom i
            ib=molstartindex+pid2 + (i * num_atoms_molecule) # index of copy atom j
            exclusions.append((ia,ib))
    return bonds
        
def storeAngles(f, types, angletypes, angletypeparams, angles, num_atoms_molecule, num_molecule_copies, molstartindex):
    line = ''
    angles_tmp = []
    pos = f.tell()
    line=f.readlastline()
    while not 'angles' in line:
        line = f.readline()
        if 'moleculetype' in line or not line:
            f.seek(pos)
            return angles
        
    line = f.readline()
    while(len(line) > 1 and not '[' in line):
        if line[0] == ";": # skip comment lines
            line = f.readline()
            continue
        tmp = line.split()
        lookup=(len(tmp)<=4)
        pid1, pid2, pid3 = map(int, tmp[0:3])
        if lookup:
            t1, t2, t3 = types[pid1-1], types[pid2-1], types[pid3-1]
            try:
                antypeid = angletypes[t1][t2][t3]
            except KeyError:
                #todo: is this good style?
                t1, t3 = t3, t1
                antypeid = angletypes[t1][t2][t3]
        else:
            #check if we need to make new type
            temptype=ParseAngleTypeParam(line) 
            antypeid=FindType(temptype, angletypeparams)
            if antypeid==None:
                antypeid=len(angletypeparams)
                angletypeparams.update({antypeid:temptype})
                
        angles_tmp.append((pid1, pid2, pid3, antypeid)) # store angletypes for this molecule
        
        line = f.readline()
        
    # extend angles to copies of this molecule
    angles_per_mol = len(angles_tmp)
    for i in range(num_molecule_copies):
        for j in range(angles_per_mol):
            pid1, pid2, pid3, antypeid = angles_tmp[j][0:4]
            ia=molstartindex+pid1 + (i * num_atoms_molecule) # index of copy atom i
            ib=molstartindex+pid2 + (i * num_atoms_molecule) # index of copy atom j
            ic=molstartindex+pid3 + (i * num_atoms_molecule) # index of copy atom k
            if antypeid in angles:
                angles[antypeid].append((ia, ib, ic))
            else:
                angles.update({antypeid:[(ia, ib, ic)]})
    return angles  

def storeDihedrals(f, types, dihedraltypes, dihedraltypeparams, dihedrals, num_atoms_molecule, num_molecule_copies, molstartindex):
    line = ''
    dihedrals_tmp = []
    pos = f.tell()
    line=f.readlastline()
    while not 'dihedrals' in line:
        line = f.readline()
        if 'moleculetype' in line or not line:
            f.seek(pos)
            return dihedrals

    line = f.readline()
    while(len(line) > 1 and not '[' in line):
        if line[0] == ";": # skip comment lines
            line = f.readline()
            continue
        tmp = line.split()
        lookup=(len(tmp)<=5)
        pid1, pid2, pid3, pid4 = map(int, tmp[0:4])
        if lookup:
            t1, t2, t3, t4 = types[pid1-1], types[pid2-1], types[pid3-1], types[pid4-1] # get types of particles
            try:
                dihtypeid = dihedraltypes[t1][t2][t3][t4]
            #if t1 not in dihedraltypes: # interactions in the other way
            except KeyError:
                t1, t2, t3, t4 = t4, t1, t2, t3
                dihtypeid = dihedraltypes[t1][t2][t3][t4]
        else:
            #check if we need to make new type
            temptype=ParseDihedralTypeParam(line)
            dihtypeid=FindType(temptype, dihedraltypeparams)
            if dihtypeid==None:
                dihtypeid=len(dihedraltypeparams)
                dihedraltypeparams.update({dihtypeid:temptype})
        
        dihedrals_tmp.append((pid1, pid2, pid3,pid4, dihtypeid)) # store angletypes for this molecule
        line = f.readline()
        
    # extend angles to copies of this molecule
    dihedrals_per_mol = len(dihedrals_tmp)
    for i in range(num_molecule_copies):
        for j in range(dihedrals_per_mol):
            pid1, pid2, pid3, pid4, dihtypeid = dihedrals_tmp[j][0:5]
            ia=molstartindex+pid1 + (i * num_atoms_molecule) # index of copy atom i
            ib=molstartindex+pid2 + (i * num_atoms_molecule) # index of copy atom j
            ic=molstartindex+pid3 + (i * num_atoms_molecule) # index of copy atom k
            id=molstartindex+pid4 + (i * num_atoms_molecule) # index of copy atom l
            if dihtypeid in dihedrals:
                dihedrals[dihtypeid].append((ia, ib, ic, id))
            else:
                dihedrals.update({dihtypeid:[(ia, ib, ic, id)]})
    return dihedrals
    

def setBondedInteractions(system, bonds, bondtypeparams, fpl=None):
    list={}
    bc=0
    for id, bondlist in bonds.iteritems():
        if (not fpl): fpl = espresso.FixedPairList(system.storage)
        fpl.addBonds(bondlist)
        bc+=len(bondlist) 
        bdinteraction=bondtypeparams[id].createEspressoInteraction(system, fpl)
        if bdinteraction:
            system.addInteraction(bdinteraction)
            list.update({id: bdinteraction})
    return list

def setAngleInteractions(system, angles, angletypeparams, fpl=None):
    list={}
    
    for id, anglelist in angles.iteritems():
        if (not fpl):fpl = espresso.FixedTripleList(system.storage)
        fpl.addTriples(anglelist)
        angleinteraction=angletypeparams[id].createEspressoInteraction(system, fpl)
        if angleinteraction:
            system.addInteraction(angleinteraction)
            list.update({id: angleinteraction})
    return list

def setDihedralInteractions(system, dihedrals, dihedraltypeparams, fpl=None):
    list={}
    
    for id, dihedrallist in dihedrals.iteritems():
        if (not fpl):fpl = espresso.FixedQuadrupleList(system.storage)
        fpl.addQuadruples(dihedrallist)
        dihedralinteraction=dihedraltypeparams[id].createEspressoInteraction(system, fpl)
        if dihedralinteraction:
            system.addInteraction(dihedralinteraction)
            list.update({id: dihedralinteraction})
    return list
        
def setLennardJonesInteractions(system, defaults, atomtypeparams, verletlist, cutoff, hadress=False, ftpl=None):
    """ Set lennard jones interactions which were read from gromacs based on the atomypes"""
    if not hadress:
	interaction = espresso.interaction.VerletListLennardJones(verletlist)
    else:
	interaction=espresso.interaction.VerletListHadressLennardJones(verletlist, ftpl)
    #interaction = espresso.interaction.VerletListLennardJonesGromacs(verletlist)
    
    print "Setting up Lennard-Jones interactions"
    if defaults:
        if int(defaults['combinationrule'])!=2:
            for atnr, at in atomtypeparams.iteritems():
                c6=float(at['sig'])
                c12=float(at['eps'])
                if c6==0: continue
                sig = pow(c12/c6,1.0/(6.0))
                eps = 0.25*c6*pow(sig,-6.0)
                at['sig']=sig
                at['eps']=eps
                print "WARNING: Converted atomtype number ", atnr, "to sigma, epsilon parameters", " sig= ", sig, " eps=", eps
    
    for i in range(len(atomtypeparams)):
        for j in range(i, len(atomtypeparams)):
            pi=atomtypeparams[i]
            pj=atomtypeparams[j]
            if pi!=pj:
                sig=0.5*(float(pi['sig'])+float(pj['sig']))
                eps=math.sqrt(float(pi['eps'])*float(pj['eps']))
            else:
                sig=float(pi['sig'])
                eps=float(pi['eps'])
            if (sig>0 and eps >0):
                #print "Setting LJ interaction for", i, j, "to sig ", sig, "eps", eps, "cutoff", cutoff
                ljpot= espresso.interaction.LennardJones(epsilon=eps, sigma=sig, shift='auto', cutoff=cutoff)
                if not hadress:
		    interaction.setPotential(type1=i, type2=j, potential=ljpot)
		else:
		    interaction.setPotentialAT(type1=i, type2=j, potential=ljpot)
    system.addInteraction(interaction)
    return interaction

def setCoulombInteractions(system, verletlist, rc, types, epsilon1, epsilon2,kappa, hadress=False, ftpl=None):
 

    pref=138.935485 # we want gromacs units, so this is 1/(4 pi eps_0) ins units of kJ mol^-1 e^-2
    
    pot = espresso.interaction.ReactionFieldGeneralized(prefactor=pref, kappa=kappa, epsilon1=epsilon1, epsilon2=epsilon2, cutoff=rc)
    #pot = espresso.interaction.CoulombTruncated(qq=0.6724, cutoff=rc, shift=0)
    if (not hadress):	
	interaction=espresso.interaction.VerletListReactionFieldGeneralized(verletlist)
    else: 
	interaction=espresso.interaction.VerletListHadressReactionFieldGeneralized(verletlist, ftpl)
   # interaction=espresso.interaction.VerletListCoulombTruncated(verletlist)
            
    for i in range(max(types)+1):
        for k in range(i, max(types)+1):
	    if (not hadress):
		interaction.setPotential(type1=i, type2=k, potential=pot)
	    else:
		interaction.setPotentialAT(type1=i, type2=k, potential=pot)

    system.addInteraction(interaction)
    return interaction

def setTabulatedInteractions(potentials, particleTypes, system, interaction):
    """Set interactions for all given particle types.
    Return value is a system with all interactions added.

    Keyword arguments:
    potentials -- is a dictionary where key is a string composed
    of two particle types and value is a potential.
    example: {"A_A":potAA, "A_B":potAB, "B_B":potBB}
    particleTypes -- is a dictionary where key is the particle type, and
    value is a list of particles of that type.
    example: {"A":["A1m", "A2m"],"B":["B1u","B2u"]}
    system -- is the system to which the interaction will be added
    interaction -- is the interaction to which to add the potentials
    """
    allparticles = []
    for k, v in particleTypes.iteritems():
        for i in v:
            allparticles.append((i,k)) # create tuples: (particle, type)
    
    for i in range(len(allparticles)):
        for j in range(i, len(allparticles)):
            type1 = allparticles[i][1]
            type2 = allparticles[j][1]
            key = type1+"_"+type2
            interaction.setPotential(i, j, potentials[key])

    system.addInteraction(interaction)
    return interaction




def convertTable(gro_in_file, esp_out_file, sigma=1.0, epsilon=1.0, c6=1.0, c12=1.0):
    """Convert GROMACS tabulated file into ESPResSo++ tabulated file (new file
    is created). First column of input file can be either distance or angle.
    For non-bonded files, c6 and c12 can be provided. Default value for sigma, epsilon,
    c6 and c12 is 1.0. Electrostatics are not taken into account (f and fd columns).
    
    Keyword arguments:
    gro_in_file -- the GROMACS tabulated file name (bonded, nonbonded, angle
    or dihedral).
    esp_out_file -- filename of the ESPResSo++ tabulated file to be written.
    sigma -- optional, depending on whether you want to convert units or not.
    epsilon -- optional, depending on whether you want to convert units or not.
    c6 -- optional
    c12 -- optional
    """

    

    # determine file type
    bonded, angle, dihedral = False, False, False
    if gro_in_file[6] == "b":
        bonded = True
    if gro_in_file[6] == "a":
        angle  = True
        bonded = True
    if gro_in_file[6] == "d":
        dihedral = True
        bonded = True

    fin = open(gro_in_file, 'r')
    fout= open(esp_out_file, 'w')

    if bonded: # bonded has 3 columns
        for line in fin:
            if line[0] == "#": # skip comment lines
                continue
            
            columns = line.split()
            r = float(columns[0])
            f = float(columns[1]) # energy
            fd= float(columns[2]) # force
            
            # convert units
            if angle or dihedral: # degrees to radians
                r = math.radians(r)
                fd=fd*180/math.pi
            else:
                r = r / sigma
            e = f / epsilon
            f = fd*sigma / epsilon
            
            if (not angle and not dihedral and r != 0) or \
                 (angle and r <= 3.1415 and r > 0) or \
                  (dihedral and r >= -3.1415 and r <= 3.1415):
                fout.write("%15.8g %15.8g %15.8g\n" % (r, e, f))
    
    else: # non-bonded has 7 columns
        for line in fin:
            if line[0] == "#": # skip comment lines
                continue
            
            columns = line.split()
            r = float(columns[0])
            #f = float(columns[1]) # electrostatics not implemented yet
            #fd= float(columns[2]) # electrostatics not implemented yet
            g = float(columns[3]) # dispersion
            gd= float(columns[4])
            h = float(columns[5]) # repulsion
            hd= float(columns[6])
            
            e = c6*g + c12*h
            f = c6*gd+ c12*hd
            
            # convert units
            r = r / sigma
            e = e / epsilon
            f = f*sigma / epsilon
            
            if r != 0: # skip 0
                fout.write("%15.8g %15.8g %15.8g\n" % (r, e, f))
    
    fin.close()
    fout.close()
