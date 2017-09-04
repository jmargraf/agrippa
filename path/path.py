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

print(mollist)
print(etot)
print(reactions)

nPaths = int(1e4)
startmol = 89
endmol = 2

# run path search
finalpaths = []

print("running path search")
for iPaths in range(nPaths):
    if iPaths%10000 == 0:
        print(iPaths/nPaths)
    Emax = -100000.0
    Estart = etot[startmol]
    molid = startmol
    iStep = 0
    #print("Path %s :" % iPaths)
    curpath = []
    lastreact = []
    available = [89,139]

    while iStep < 20:
        molid = random.choice(available)
        decomp = []
        contingent = []
        flagdeco = False
        flagcont = False

        # search possible reactions
        #print("    search possible reactions")
        for react in reactions:
            if react[0] == molid:
                decomp.append(react)
            elif (react[1] == molid) or (react[2] == molid):
                contingent.append(react)

        # select reaction
        #print("    select reaction")
        reactid=random.randrange(len(decomp)+len(contingent))
        if reactid > (len(decomp)-1):
            react = contingent[reactid-len(decomp)]
            flagcont = True
        else:
            react = decomp[reactid]
            flagdeco = True


        # check reactions
        #print("    check reactions")
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
        else:
            lastreact = react

        # add to path and select next molid
        #print("    add to path")
        if flagdeco:
            dG = (-etot[react[0]] + etot[react[1]] + etot[react[2]])#*627.5
            if dG > Emax:
                Emax = dG 

            reastring =("step " + str(iStep) + ": " + 
                        mollist[react[0]] + " (" + str(react[0]) + ") ->  " +
                        mollist[react[1]] + " (" + str(react[1]) + ") + " + 
                        mollist[react[2]] + " (" + str(react[2]) + ")  | dG =  " + str(dG) )
            #print(reastring)
            curpath.append([reastring,dG,"deco",react[0],react[1],react[2]])            
            available.append(react[1])
            available.append(react[2])
     
            if (react[1] == endmol) or (react[2] == endmol):
                #print("    ! path %s reached target, Emax = %s kcal mol-1" % (iPaths,Emax))
                finalpaths.append([curpath,Emax])
                break
            else:
                #if (react[1]==20) or (react[1]==21):
                #    molid = react[2]
                #elif (react[2]==20) or (react[2]==21):
                #    molid = react[1]
                #else:
                molid = react[random.randint(1,2)]        
            
        elif flagcont:
            dG = (+etot[react[0]] - etot[react[1]] - etot[react[2]])#*627.5
            if dG > Emax:
                Emax = dG

            reastring =("step " + str(iStep) + ": " + 
                        mollist[react[1]] + " (" + str(react[1]) + ") +  " +
                        mollist[react[2]] + " (" + str(react[2]) + ") -> " +
                        mollist[react[0]] + " (" + str(react[0]) + ")  | dG =  " + str(dG) )
            #print(reastring)
            curpath.append([reastring,dG,"cont",react[1],react[2],react[0]])       
            available.append(react[0])

            if (react[0] == endmol):
                #print("    ! path %s reached target, Emax = %s kcal mol-1" % (iPaths, Emax))
                finalpaths.append([curpath,Emax])
                break
            else:
                molid = react[0]

        iStep += 1

# print examples of best paths
sortedpaths = sorted(finalpaths, key=lambda path: path[1])

lastbar = 1000000.0
printtop = 5
iprint = 0
print("printing low barrier example paths:")
for path in sortedpaths:
    if iprint >= printtop:
        break
    elif path[1] == lastbar:
        continue
    else:
        lastbar = path[1]
        for reaction in path[0]:
            print(reaction[0])
        print("    ! Emax = %s eV" % (path[1]))
        iprint += 1


print("printing all lowest barrier reactions")
print("    ! Emax = %s eV" % (sortedpaths[0][1]))
lowbar = []
for path in sortedpaths:
    if path[1] > sortedpaths[0][1]:
        break
    else:
        for reaction in path[0]:
            oldflag = False
            for oldreaction in lowbar:
                if reaction[1:5] == oldreaction[1:5]:
                    oldflag = True
            if not oldflag:
                print(reaction[0])
                print(reaction[2],reaction[3],reaction[4],reaction[5])
                lowbar.append(reaction)

print("printing shortest path for lowest barrier reaction:")
print("    ! Emax = %s eV" % (sortedpaths[0][1]))
minlen = 100
for path in sortedpaths:
    if path[1] > sortedpaths[0][1]:
        break
    else:
        if len(path[0])<minlen:
            shortest = path[0]
            minlen = len(path[0])

for reaction in shortest:
    print(reaction[0])


