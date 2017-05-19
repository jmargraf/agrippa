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
            etot.append(float(dat[1]))

    reactions = []

    # read reactions 
    with open('react.txt','r') as f:
        for line in f:          
            dat = line.split()
            reactions.append([int(dat[2]),int(dat[5]),int(dat[8])])

    return mollist,etot,reactions


print("reading reaction network")
mollist,etot,reactions = Read_Reactions()

#print(mollist)
#print(etot)
#print(reactions)

nPaths = 1000000
startmol = 24
endmol = 2

# run path search
finalpaths = []

print("running path search")
for iPaths in range(nPaths):
    Emax = -100000.0
    Estart = etot[startmol]
    molid = startmol
    iStep = 0
    #print("Path %s :" % iPaths)
    curpath = []
    lastreact = []
    available = [20,21]
    
    while iStep < 10:
        decomp = []
        contingent = []
        flagdeco = False
        flagcont = False

        # search possible reactions
        for react in reactions:
            if react[0] == molid:
                decomp.append(react)
            elif (react[1] == molid) or (react[2] == molid):
                contingent.append(react)

        # select reaction
        reactid=random.randrange(len(decomp)+len(contingent))
        if reactid > (len(decomp)-1):
            react = contingent[reactid-len(decomp)]
            flagcont = True
        else:
            react = decomp[reactid]
            flagdeco = True


        # check reactions
        if(react == lastreact):
            continue
        elif(flagcont):
            avail = False
            for availmol in available:
                if molid == react[1]:
                    if react[2] == availmol:
                        avail = True
                elif molid == react[2]:
                    if react[1] == availmol:
                        avail = True
            if not avail:
                continue
        else:
            lastreact = react

        # add to path and select next molid
        if flagdeco:
            dG = (-etot[react[0]] + etot[react[1]] + etot[react[2]])*627.5
            if dG > Emax:
                Emax = dG 

            reastring =("step " + str(iStep) + ": " + 
                        mollist[react[0]] + " (" + str(react[0]) + ") ->  " +
                        mollist[react[1]] + " (" + str(react[1]) + ") + " + 
                        mollist[react[2]] + " (" + str(react[2]) + ")  | dG =  " + str(dG) )
            #print(reastring)
            curpath.append(reastring)            
            available.append(react[1])
            available.append(react[2])
     
            if (react[1] == endmol) or (react[2] == endmol):
                #print("    ! path %s reached target, Emax = %s kcal mol-1" % (iPaths,Emax))
                finalpaths.append([curpath,Emax])
                break
            else:
                if (react[1]==20) or (react[1]==21):
                    molid = react[2]
                elif (react[2]==20) or (react[2]==21):
                    molid = react[1]
                else:
                    molid = react[random.randint(1,2)]        
            
        elif flagcont:
            dG = (+etot[react[0]] - etot[react[1]] - etot[react[2]])*627.5
            if dG > Emax:
                Emax = dG

            reastring =("step " + str(iStep) + ": " + 
                        mollist[react[1]] + " (" + str(react[1]) + ") +  " +
                        mollist[react[2]] + " (" + str(react[2]) + ") -> " +
                        mollist[react[0]] + " (" + str(react[0]) + ")  | dG =  " + str(dG) )
            #print(reastring)
            curpath.append(reastring)
            available.append(react[0])

            if (react[0] == endmol):
                #print("    ! path %s reached target, Emax = %s kcal mol-1" % (iPaths, Emax))
                finalpaths.append([curpath,Emax])
                break
            else:
                molid = react[0]

        iStep += 1

# print best paths
sortedpaths = sorted(finalpaths, key=lambda path: path[1])

lastbar = 1000000.0
printtop = 10
iprint = 0
print("printing lowest barrier paths:")
for path in sortedpaths:
    if iprint >= printtop:
        break
    elif path[1] == lastbar:
        continue
    else:
        lastbar = path[1]
        for reaction in path[0]:
            print(reaction)
        print("    ! Emax = %s kcal mol-1" % (path[1]))
        iprint += 1








