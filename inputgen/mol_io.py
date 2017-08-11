#!/usr/bin/env python

from __future__ import print_function


#       xyz,atom,charge,spin = io.ReadMol(FileName)
#       xyz,atom,charge,spin = io.ReadXYZ(FileName)
#       io.OrcaIn(keywords,atom,xyz,charge,spin)


###################################################
# Subroutines
# generate backbones (linear and branched)

def ReadMol(FileName):
    name = FileName + ".mol"
    nAtoms = 0
    xyz = []
    atom = []
    charge = 0
    spin = 1
    with open(name) as fi:
        f = []
        for line in fi:
            f.append(line)
        data = f[3].strip().split()
        nAtoms = int(data[0])
        for i in range(nAtoms):
            data = f[4+i].strip().split()
            atom.append(data[3])
            xyz.append([float(data[0]),float(data[1]),float(data[2])])
        for count in range(len(f)):
           data = f[count].strip().split()
           if len(data)>1 and data[1] == "RAD":
               nrad=int(data[2])
               for i in range(nrad):
                   spin += (int(data[4+2*i])-1)
           elif len(data)>1 and data[1] == "CHG":
               nchg=int(data[2])
               for i in range(nchg):
                   charge += int(data[4+2*i])

    return xyz,atom,charge,spin

def ReadXYZ(FileName):
#   print("placeholder for read xyz")
    name = FileName + ".xyz"
    nAtoms = 0
    xyz = []
    atom = []
    charge = 0
    spin = 1
    with open(name) as fi:
        f = []
        for line in fi:
            f.append(line)
        data = f[0].strip().split()
        nAtoms = int(data[0])
        data = f[1].strip().split()
        if len(data)==2:
            charge = int(data[0])
            spin = int(data[1])
        for i in range(nAtoms):
            data = f[2+i].strip().split()
            atom.append(data[0])
            xyz.append([float(data[1]),float(data[2]),float(data[3])])

    return xyz,atom,charge,spin

def ReadAimsGeo():
    name = 'geometry.in'
    #nAtoms = 0
    atoms = []
    xyz = []
    with open(name) as fi:
        f = []
        for line in fi:
            tmp = line.strip()
            if tmp == "":
                break
            elif tmp[0] == "#":
                continue
            else:
                f.append(tmp)
        for line in f:
            atom = line.split()
            xyz.append([float(atom[1]),float(atom[2]),float(atom[3])])
            atoms.append(atom[4])
    return xyz,atoms


def ReadXYZtraj(FileName):
    name = FileName + ".xyz"
    nAtoms = 0
    xyz = []
    atom = []
    with open(name) as fi:
        f = []
        for line in fi:
            f.append(line)
    data = f[0].strip().split()
    nAtoms = int(data[0])
    nFrames = len(f)/(nAtoms+2)
#   print(nFrames)
    for iFrame in range(nFrames):
        for i in range(nAtoms):
            data = f[i+2+iFrame*(nAtoms+2)].strip().split()
            #print(data)
            atom.append(data[0])
            xyz.append([float(data[1]),float(data[2]),float(data[3])])

    return xyz,atom,nAtoms



def OrcaIn(FileName,keywords,atom,xyz,charge,spin):
#   print("placeholder for write orca")
    name =  FileName + ".inp"
    print("# orca input made with inputgen.py",file=open(name,'w'))

    for text in keywords:
        print(text,file=open(name,'a'))

    print("\n",file=open(name,'a'))
    text = "* xyz " + str(charge) + " " + str(spin)
    print(text,file=open(name,'a'))
    for i in range(len(atom)):
        text = atom[i] + "  " + str(xyz[i][0]) + "  " + str(xyz[i][1]) + "  " + str(xyz[i][2]) 
        print(text,file=open(name,'a'))
    print("*",file=open(name,'a'))
    print("\n",file=open(name,'a'))


def PrintXYZ(name,atom,xyz):
    nAtoms = len(atom)
    print(nAtoms)
    print(name)
    for i in range(nAtoms):
        text = atom[i] + "  " + str(xyz[i][0]) + "  " + str(xyz[i][1]) + "  " + str(xyz[i][2])
        print(text)


def PrintAimsGeo(name,atom,xyz):
    nAtoms = len(atom)
    text = "# " + name
    print(text)
    for i in range(nAtoms):
        text = "    atom    " + str(xyz[i][0]) + "  " + str(xyz[i][1]) + "  " + str(xyz[i][2]) + "  " + atom[i] 
        print(text)


def MergeXYZ(atom1,xyz1,atom2,xyz2):
    
    atom = []
    xyz = []
    frags = []

    for coord in xyz1:
        xyz.append(coord)
        frags.append(1)
    for coord in xyz2:
        xyz.append(coord)
        frags.append(2)
    for label in atom1:
        atom.append(label)
    for label in atom2:
        atom.append(label)


    return atom,xyz,frags
