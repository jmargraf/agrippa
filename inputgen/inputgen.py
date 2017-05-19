#!/usr/bin/env python

from __future__ import print_function
import sys, getopt
import mol_io as io
import react

def main(argv):

    # Initial Input
    FileName = " " 
    pro1Name = " "
    pro2Name = " "
    GeoType = "mol"
    Template = "template.txt"
    Code = "orca"
    Reaction = False

    xyz = []
    atom = []
    charge = 0
    spin = 1

#   python inputgen.py -i [input] --mol --orca --temp [template.txt]

    try:
        opts,args = getopt.getopt(argv,"hi:",['mol','xyz','temp=','orca','pro1=','pro2=','react'])
    except getopt.GetoptError:
        print("inputgen.py -i [input] --mol --orca --temp [template.txt]")
    for opt, arg in opts:
        if opt == "-h":
            print("inputgen.py -i [input] {--mol --orca --temp [template.txt]}")
            print("    -i           : input name (w/o .mol/.xyz)")
            print("    --mol        : read .mol file (default)")
            print("    --xyz        : read .xyz file")
            print("    --orca       : write orca input file")
            print("    --temp       : read from specific template file (default is template.txt)")
            print("    --react      : generate input for reaction")
            print("    --pro1       : name of first product")
            print("    --pro2       : name of second product")
            print("    -h           : print this help")
            print("")
            sys.exit()
        elif opt == "-i":
            FileName = arg
        elif opt == "--mol":
            GeoType = "mol"
        elif opt == "--temp":
            Template = arg
        elif opt == "--xyz":
            GeoType = "xyz"
        elif opt == "--react":
            Reaction = True
        elif opt == "--pro1":
            pro1Name = arg
        elif opt == "--pro2":
            pro2Name = arg
 
    if not Reaction: # Standard input from single geo

        # open input file
        if GeoType == "mol":
            xyz,atom,charge,spin = io.ReadMol(FileName)
        elif GeoType == "xyz":
            xyz,atom,charge,spin = io.ReadXYZ(FileName)

    if Reaction: # Generate reaction geometry first

        xyz,atom,charge,spin = react.GenReaction(FileName,pro1Name,pro2Name)

        exit()

    print(charge)
    print(spin)
    print(atom)
    print(xyz)

    # read template
    with open(Template) as temp:
        keywords = temp.read().splitlines()

    print(keywords)

    # write input
    if Code == "orca":
        io.OrcaIn(FileName,keywords,atom,xyz,charge,spin)

if __name__ == "__main__":
   main(sys.argv[1:])


