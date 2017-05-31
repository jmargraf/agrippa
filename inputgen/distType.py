#!/usr/bin/env python

from __future__ import print_function
import sys
import mol_io as io
import numpy as np

def Set_Bondlength(xyz,iA,iB,shift):



if len(sys.argv) != 4:
    print("Provide filename and atomnumbers")
    quit()

name = sys.argv[1]
atomA = sys.argv[2]
atomB = sys.argv[3]

xyz,atom,nAtoms = io.ReadXYZtraj(name)

nFrame = len(xyz)/nAtoms

for iFrame in range(nFrame):

    for iA in range(nAtoms):
        if atom[iA+iFrame*nAtoms] != atomA:
            continue
        for iB in range(nAtoms):
            if atom[iB+iFrame*nAtoms] != atomB:
                continue            
            elif iB == iA:
                continue
            coordA = np.array(xyz[iA+iFrame*nAtoms])
            coordB = np.array(xyz[iB+iFrame*nAtoms])

            vector = coordA-coordB

            distance = np.linalg.norm(vector)

            print(distance)

