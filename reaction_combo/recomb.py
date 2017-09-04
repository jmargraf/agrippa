#! /usr/bin/env python

from __future__ import division, absolute_import, print_function
import numpy as np
import math
import random

def Read_Reactions():
    mollist = []
    etot = []

    # read molecules and energies
    with open('mols.txt','r') as f:
        for line in f:    
            dat = line.split()
            mollist.append(dat[0])
            etot.append(0.0)#float(dat[1]))

    reactions = []

    # read reactions 
    with open('react.txt','r') as f:
        for line in f:          
            dat = line.split()
            reactions.append([int(dat[2]),141,int(dat[5]),int(dat[8]),dat[9]])

    return mollist,etot,reactions


def Duplicate(re,refre):
    dupe = False
    for ref in refre:
        if ((sorted(re[0:2])==sorted(ref[0:2]))and
              (sorted(re[2:4])==sorted(ref[2:4]))):
            dupe = True
            break
        elif ((sorted(re[0:2])==sorted(ref[2:4]))and
              (sorted(re[2:4])==sorted(ref[0:2]))):
            dupe = True
            break

    return dupe


#print("reading reaction network")
mollist,etot,reactions = Read_Reactions()

#print(mollist)
#print(etot)
#print(reactions)

count = len(reactions)-1

new_reactions = []

# search for transfer reactions
for react1 in reactions:
    for react2 in reactions:
        re = None
        if react1[0] == react2[0]:
            # possibly a substitution?

            continue
        # rearrangement
        if sorted(react1[2:4]) == sorted(react2[2:4]):
            re = [react1[0],141,react2[0],141]
            rclass = 'rea:'+react1[4]+react2[4]

        #transfer
        elif sorted(react1[2:4])[0] == sorted(react2[2:4])[0]:
            re = [react1[0],sorted(react2[2:4])[1],react2[0],sorted(react1[2:4])[1]]
            rclass = 'tra:'+react1[4]+react2[4]
        elif sorted(react1[2:4])[1] == sorted(react2[2:4])[1]:
            re = [react1[0],sorted(react2[2:4])[0],react2[0],sorted(react1[2:4])[0]]
            rclass = 'tra:'+react1[4]+react2[4]


        if re != None:
            if ((not Duplicate(re,reactions)) and (not Duplicate(re,new_reactions))):
                new_reactions.append(re)
                count += 1
                print(count,mollist[re[0]]," + ",mollist[re[1]]," -> ",mollist[re[2]]," + ",mollist[re[3]],rclass)


