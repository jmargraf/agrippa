#!/usr/bin/env python

from __future__ import print_function
import sys
import mol_io as io
import numpy as np

if len(sys.argv) != 4:
    print("Provide filename and atomnumbers")
    quit()

name = sys.argv[1]
atomA = int(sys.argv[2])
atomB = int(sys.argv[3])

xyz,atom,nAtoms = io.ReadXYZtraj(name)

nFrame = len(xyz)/nAtoms

for iFrame in range(nFrame):

    coordA = np.array(xyz[atomA-1+iFrame*nAtoms])
    coordB = np.array(xyz[atomB-1+iFrame*nAtoms])

    vector = coordA-coordB

    distance = np.linalg.norm(vector)

    print(distance)

