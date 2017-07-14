import numpy as np
from scipy.optimize import least_squares
import mol_io as io

def read_ref(refname):
    with open(refname,'r') as f:
        data = list(f)
    mol = np.empty(len(data))
    ref = np.empty(len(data))
    for i,line in enumerate(data):
        tmp = line.split()
        mol[i] = float(i) 
        ref[i] = float(tmp[2])
    return mol,ref

def runrefcalcs(x):
    Eref = {}
    with open(refmols,'r') as f:
        data = list(f)
    for line in data:
        refmol = line.split()
        #print(refmol[0])
        Eref[refmol[0]] = runorca(refmol[0],1,x)
    return Eref

def runmolcalcs(x,mol,Eref):
    result = np.zeros(mol.size)    
    with open(refname,'r') as f:
        data = list(f)
    for i,line in enumerate(data):
        calc = line.split()
        name = calc[0]
        comp = calc[3]
        mult = int(calc[4])
        corr = ErelCorrect(comp,Eref)
        result[i] = runorca(name,mult,x)-corr
        #print(name,result[i])
    return result

def ErelCorrect(comp,Eref):
    A = 0.0
    muH = Eref['h2']/2.0
    for letter in comp:
        if letter == 'h':
            A += muH
        elif letter == 'o':
            A += Eref['h2o'] - 2.0*muH 
        elif letter == 'c':
            A += Eref['co'] - (Eref['h2o'] - 2.0*muH)
        elif letter == 'n':
            A += Eref['nh3'] - 3.0*muH
        elif letter == 'b':
            A += Eref['bh3'] - 3.0*muH
        elif letter == 'f':
            A += Eref['f2']/2.0
        elif letter == 'a':  #Warning, actuall Al
            A += Eref['alh3'] - 3.0*muH
        elif letter == 'l':  # Cl
            A += Eref['cl2']/2.0
        elif letter == 'i':  # Si! This is obviously not ideal...
            A += Eref['sih4'] - 4.0*muH
        elif letter == 's':
            A += Eref['h2s'] - 2.0*muH
        elif letter == 'p':
            A += Eref['ph3'] - 3.0*muH
    return A

def runorca(name,mult,par):
    import subprocess
    from subprocess import Popen
    Etot = 0.0
    xyz,atom,charge,spin = io.ReadXYZ(name)
    keywords  = '! PBE MINIX def2/J \n'
    parstring = ("%method \n"+
                 "  XKappa  " + str(par[0]) + "\n" +
                 "  XMuePBE " + str(par[1]) + "\n" +
                 "end \n")
    io.OrcaIn(name,keywords,parstring,atom,xyz,0,mult)
    inpname = name + '.inp'
    outname = name + '.out'

    proc = Popen(['orca',inpname],stdout=subprocess.PIPE)

    for line in proc.stdout.readlines():
        data=line.split()
        if len(data)>2:
            if (data[0] == 'FINAL') and (data[1] == 'SINGLE'):
                Etot = float(data[4])

    #with open(filename,'r') as f:
    #    data = list(f)
    #tmp = data[0].split()
    #Etot = par[0]*float(tmp[0]) + par[1]

    return Etot

def model(x,mol):
# calculate datapoint
    result = np.zeros(mol.size)    
    #for o,par in enumerate(x):
    #   result += par*mol**o  
    Eref = runrefcalcs(x)
    #print('Eref',Eref)
    result = runmolcalcs(x,mol,Eref) 
    return result

def fun(x,mol,ref):
# calculate error
    result = (model(x,mol))*627.5
    #print(result,ref)
    return result - ref

refname = 'molrefs.txt'
refmols = 'energyrefs.txt'
mol, ref = read_ref(refname)
x0 = np.array([0.80,0.20000])

#print(mol)
#print(ref)

res = least_squares(fun, x0, jac='2-point', bounds=(0,10), args=(mol,ref), verbose=2,diff_step=0.001)

print("final parameters: ",res.x)
#print("final residuals: ",res.fun)
molnames = []
with open(refname,'r') as f:
    data = list(f)
    for i,line in enumerate(data):
        tmp = line.split()
        molnames.append(tmp[0])

for i,residual in enumerate(res.fun):
    print(molnames[i],residual)

#error = fun(res.x,mol,ref)

