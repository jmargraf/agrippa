#!/usr/bin/env python

from __future__ import print_function
from rdkit import Chem
from rdkit.Chem import RWMol
from rdkit.Chem import AllChem
import mol_io as io
import topomod
import math
import random

# xyz,atom,charge,spin = react.GenReaction(FileName,pro1Name,pro2Name)

def GenReaction(EdName,Pro1Name,Pro2Name):
    EdFile = EdName + ".mol"
    Pro1File = Pro1Name + ".mol"
    Pro2File = Pro2Name + ".mol"

    edmol   = Chem.MolFromMolFile(EdFile,sanitize=True,removeHs=False)
    pro1mol = Chem.MolFromMolFile(Pro1File,sanitize=True,removeHs=False)
    pro2mol = Chem.MolFromMolFile(Pro2File,sanitize=True,removeHs=False)
    
    print(EdFile)
    print(Pro1File)
    print(Pro2File)

#   print(Chem.MolToMolBlock(edmol))
#   print(Chem.MolToMolBlock(pro1mol))
#   print(Chem.MolToMolBlock(pro2mol))

    xyzed,atomed,charge,spin = io.ReadMol(EdName)
    xyzpro1,atompro1,charge,spin = io.ReadMol(Pro1Name)
    xyzpro2,atompro2,charge,spin = io.ReadMol(Pro2Name)

    xyzed_new = []
    atomed_new = []

#   io.PrintXYZ('Educt',atomed,xyzed)

    # Align first product to molecule
    bestpair = []
    bestrmsd = 100
    ifit = 0

    while ifit < 10000: # repeatedly tests random atom-pairlists. probably not the smartest way
        ifit += 1
        pairlist1 = range(len(atompro1))
        atomlist = range(len(atompro1))
        random.shuffle(atomlist)
        flags = [False]*len(atomed)
        weightlist = [1]*len(atompro1)
#       print(atomlist)

        for j in atomlist:
            for i in range(len(atomed)):
                if(atomed[i] != atompro1[j]):
                   continue
                elif not flags[i]:
                    pairlist1[j] = [j,i]
                    flags[i] = True
                    break

#       print(pairlist1)

        if len(pairlist1) == 1:
            # if the fragment is a single atom rd-kits AlignMol fails. Simply move atom to corresponding location
            xyzpro1[0] = xyzed[pairlist1[0][1]]
            xyzed_new.append(xyzed[pairlist1[0][1]])
            atomed_new.append(atomed[pairlist1[0][1]])
            ifit = 1000000
        else:
            for i in range(len(pairlist1)):
                weightlist[i] = pro1mol.GetAtomWithIdx(i).GetAtomicNum()**2
            rmsd = AllChem.AlignMol(pro1mol,edmol,atomMap=pairlist1,weights=weightlist)
#           print(rmsd)
            if rmsd<bestrmsd:
                bestrmsd = rmsd
                bestpair = pairlist1

    if len(pairlist1) > 1:
        # realign with best pairlist. also, write to new educt geo
        pairlist1 = bestpair
        for i in range(len(pairlist1)):
            weightlist[i] = pro1mol.GetAtomWithIdx(i).GetAtomicNum()**2
            xyzed_new.append(xyzed[pairlist1[i][1]])
            atomed_new.append(atomed[pairlist1[i][1]])
        rmsd = AllChem.AlignMol(pro1mol,edmol,atomMap=pairlist1,weights=weightlist)
        print(rmsd)
        print(Chem.MolToMolBlock(pro1mol),file=open('pro1mol.mol','w+'))
        xyzpro1,atompro1,charge,spin = io.ReadMol('pro1mol')

#   io.PrintXYZ('Product1',atompro1,xyzpro1)

#   print(Chem.MolToMolBlock(edmol),file=open('edmol.mol','w+'))
#   print(Chem.MolToMolBlock(pro1mol),file=open('pro1mol.mol','w+'))


    # Align second product to remaining framework. Same as above
    bestpair = []
    bestrmsd = 100
    ifit = 0

    while ifit < 10000:
        ifit += 1
        pairlist2 = range(len(atompro2))
        atomlist = range(len(atompro2))
        random.shuffle(atomlist)
        flags = [False]*len(atomed)
        weightlist = [1]*len(atompro2)
#       print(atomlist)

        for j in atomlist:
            for i in range(len(atomed)):
                if(atomed[i] != atompro2[j]):
                   continue
                elif any(x[1] == i for x in pairlist1): # skip atoms already used in prev step.
                    continue
                elif not flags[i]:
                    pairlist2[j] = [j,i]
                    flags[i] = True
                    break

#       print(pairlist1)

        if len(pairlist2) == 1:
            xyzpro2[0] = xyzed[pairlist2[0][1]]
            xyzed_new.append(xyzed[pairlist2[0][1]])
            atomed_new.append(atomed[pairlist2[0][1]])
            ifit = 1000000
        else:
            for i in range(len(pairlist2)):
                weightlist[i] = pro2mol.GetAtomWithIdx(i).GetAtomicNum()**2
            rmsd = AllChem.AlignMol(pro2mol,edmol,atomMap=pairlist2,weights=weightlist)
#           print(rmsd)
            if rmsd<bestrmsd:
                bestrmsd = rmsd
                bestpair = pairlist2

    if len(pairlist2) > 1:
        pairlist2 = bestpair
        for i in range(len(pairlist2)):
            weightlist[i] = pro2mol.GetAtomWithIdx(i).GetAtomicNum()**2
            xyzed_new.append(xyzed[pairlist2[i][1]])
            atomed_new.append(atomed[pairlist2[i][1]])
        rmsd = AllChem.AlignMol(pro2mol,edmol,atomMap=pairlist2,weights=weightlist)
        print(rmsd)
        print(Chem.MolToMolBlock(pro2mol),file=open('pro2mol.mol','w+'))
        xyzpro2,atompro2,charge,spin = io.ReadMol('pro2mol')

#   io.PrintXYZ('Product2',atompro2,xyzpro2)

    # shift along vector connecting fragments' centers of mass
    xyzpro1s,xyzpro2s = CoM_shift(atompro1,xyzpro1,atompro2,xyzpro2,4.0)

    # merge aligned and shifted product geometries
    atom_diss,xyz_diss,frags = io.MergeXYZ(atompro1,xyzpro1,atompro2,xyzpro2)
    atom_shift,xyz_shift,frags = io.MergeXYZ(atompro1,xyzpro1s,atompro2,xyzpro2s)

    io.PrintXYZ('Reordered Educt',atomed_new,xyzed_new)
#   io.PrintXYZ('Merged Products',atom_diss,xyz_diss)
    io.PrintXYZ('Shifted Products',atom_shift,xyz_shift)

    io.PrintAimsGeo('Reordered Educt',atomed_new,xyzed_new)
    io.PrintAimsGeo('Shifted Products',atom_shift,xyz_shift)



    return xyzed,atomed,charge,spin


def CoM_shift(atom1,xyz1,atom2,xyz2,shift):
    #calculate center of mass 1
    sumM = 0.0
    CoM1 = [0.0, 0.0, 0.0]

    for i in range(len(xyz1)):
        mi = Mass(atom1[i])
        sumM += mi
        CoM1[0] += xyz1[i][0]*mi
        CoM1[1] += xyz1[i][1]*mi
        CoM1[2] += xyz1[i][2]*mi

    CoM1[0] /= sumM
    CoM1[1] /= sumM         
    CoM1[2] /= sumM         

    #calculate center of mass 1
    sumM = 0.0
    CoM2 = [0.0, 0.0, 0.0]

    for i in range(len(xyz2)):
        mi = Mass(atom2[i])
        sumM += mi
        CoM2[0] += xyz2[i][0]*mi
        CoM2[1] += xyz2[i][1]*mi
        CoM2[2] += xyz2[i][2]*mi

    CoM2[0] /= sumM    
    CoM2[1] /= sumM         
    CoM2[2] /= sumM         

    # normalized vector from 1 to 2
    vect = [CoM1[0]-CoM2[0],CoM1[1]-CoM2[1],CoM1[2]-CoM2[2]]
    norm = math.sqrt(vect[0]**2+vect[1]**2+vect[2]**2)
    vect[0] /= norm
    vect[1] /= norm
    vect[2] /= norm
    norm = math.sqrt(vect[0]**2+vect[1]**2+vect[2]**2)
    print(norm)

    # calculate shifted coordinates
    xyz1_new = []
    for coords in xyz1:
#       xyz1_new.append([coords[0] + vect[0]*shift,coords[1] + vect[1]*shift,coords[2] + vect[2]*shift])
        xyz1_new.append([coords[0],coords[1],coords[2]])
    xyz2_new = []
    for coords in xyz2:
        xyz2_new.append([coords[0] - vect[0]*shift,coords[1] - vect[1]*shift,coords[2] - vect[2]*shift])

    return xyz1_new,xyz2_new



def Mass(El):
    if(El=='H'):
        return 1.01
    elif(El=='He'):
        return 4.0
    elif(El=='Li'):
        return 6.94
    elif(El=='Be'):
        return 9.01
    elif(El=='B'):
        return 10.81
    elif(El=='C'):
        return 12.01
    elif(El=='N'):
        return 14.01
    elif(El=='O'):
        return 16.00
    elif(El=='F'):
        return 19.00
    elif(El=='Ne'):
        return 20.18
    else:
        print("mass unknown")




