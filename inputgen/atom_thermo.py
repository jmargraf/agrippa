#!/usr/bin/env python

from __future__ import print_function
import sys
import math 
import numpy as np

if len(sys.argv) != 5:
    print("Provide filename and atomnumbers")
    quit()

Temp = float(sys.argv[1])/3.1577464e5
Pres = float(sys.argv[2])*101325.0/2.9421912e13
Mass = float(sys.argv[3])*1822.888486192 # *1822.8603451550191 #1822.888486192
Mult = float(sys.argv[4])

kb = 3.166811429E-6
kb = 1.0
h = 2.0*math.pi

Strans = kb*(math.log((((2.0*math.pi*Mass*kb*Temp)/(h*h))**(3.0/2.0))*((kb*Temp)/Pres)) + 5.0/2.0) * Temp

Sel = kb*(math.log(Mult)) * Temp

Utrans = 3.0/2.0 * kb * Temp
Htherm = kb * Temp

AtomCorr = Strans+Sel+Utrans+Htherm

print("S_trans*T = %s Ha" % Strans)
print("S_el*T    = %s Ha" % Sel)
print("U_trans   = %s Ha" % Utrans)
print("H_therm   = %s Ha" % Htherm)

print("Tot. corr = %s Ha" % AtomCorr)


