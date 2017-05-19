#! /usr/bin/env python

from __future__ import division, absolute_import, print_function
import numpy as np
import math

def MolZ(AtomType):

    if AtomType == "H":
        return 1.0
    elif AtomType == "C":
        return 6.0
    elif AtomType == "N":
        return 7.0
    elif AtomType == "O":
        return 8.0
    elif AtomType == "S":
        return 16.0
    else:
        print("unknown element %s",AtomType)

def Read_CM():
    mollist = []
    ae = []

    # read and generate all coulomb matrices
    with open('dsgdb7ae2.xyz','r') as f:
        for imol in range(7102):    
            nAtoms = int(f.readline())
            dat = f.readline().split()
            ae.append(float(dat[1]))
            atoms = []
            xyz = np.zeros((nAtoms,3))
            for i in range(nAtoms): 
                dat = f.readline().split()
                atoms.append(dat[0])
                xyz[i][0] = float(dat[1])
                xyz[i][1] = float(dat[2])
                xyz[i][2] = float(dat[3])
            # calculate coulomb matrix
            M = np.zeros((23,23))
            for i in range(nAtoms):
                for j in range(nAtoms):
                    if i == j:
                        M[i,j]=0.5*MolZ(atoms[i])**2.4
                    else:
                        M[i,j]=MolZ(atoms[i])*MolZ(atoms[j])/np.linalg.norm(xyz[i]-xyz[j])
            # sort coulomb matrix
            indexlist = np.argsort(-np.linalg.norm(M,axis=1))
            sym_sortedM=M[indexlist][:,indexlist]

            # store lower triangle as vector
            molvect = np.zeros((23*23-23)/2+23)
            iat=0
            for i in range(23):
                for j in range(i+1):
                    molvect[iat] = sym_sortedM[i][j] 
                    iat += 1

            mollist.append(molvect)

#    print(mollist[0:10])

    return mollist,ae

def KroD(i,j):
    if i==j:
        return 1.0
    else:
        return 0.0
        
def Kernel(x,y,sigma,Type):
    if Type == 1:
        norm2 = np.linalg.norm(x-y)*np.linalg.norm(x-y)
        kxy = math.exp(-norm2/(2.0*sigma*sigma))
        return kxy
    elif Type == 2:
        norm1 = np.linalg.norm(x-y,ord=1)
        kxy = math.exp(-norm1/sigma)
        return kxy

def ForwardBackwardSub(U,y):
    alpha = np.zeros(len(y))
    n = len(y)
    for i in range(n):
        v = y[i]
        for j in range(i):
            v -= U[j][i]*alpha[j]
        alpha[i] = v/U[i][i]

    for i in range(n-1,-1,-1):
        v = alpha[i]
        for j in range(n-1,i,-1):
            v -= U[i][j]*alpha[j]
        alpha[i] = v/U[i][i]
    return alpha


def MAD(est,ref):
    mae = 0.0
    print(len(est))
    print(len(ref))
    for i in range(len(est)):
        mae += abs(ref[i] - est[i]) 
    mae /= float(len(est))
    return mae    

CoulombMats = []
AEs = []

print("reading coulmob matrices")
CoulombMats, AEs = Read_CM()

print(len(CoulombMats))
print(len(AEs))

training = []
prediction = []
trainAEs = []
predictAEs = []

k=59
kmax=7102

# define training set and prediction set
for i in range(k):
    training.append(CoulombMats[i])
    trainAEs.append(AEs[i])

for i in range(k,kmax):
    if i%7 == 0: 
        training.append(CoulombMats[i]) 
        trainAEs.append(AEs[i])
    else:
        prediction.append(CoulombMats[i])
        predictAEs.append(AEs[i])


#print(len(training))
#print(len(prediction))
#print(len(training[0]))

# build kernel matrix for training
K = np.zeros((len(training),len(training)))

#Hyperparameters:
Lambda = 10.0**(-8.0) # 10.0**(-6.5)
Sigma = 2.0**12.0 #724.0
KernelType = 2

print("building K")

for i in range(len(training)):
    for j in range(len(training)):
        K[i][j] = Kernel(training[i],training[j],Sigma,KernelType) + Lambda*KroD(i,j)

#print(len(K))

# Cholesky decomp of kernel matrix 
#U = np.zeros((len(training),len(training)))
print("Cholesky decomp of K")
U = np.linalg.cholesky(K)

#print(len(U))

#alpha = ForwardSub(np.transpose(U),trainAEs)
#alpha = BackSub(U,alpha)
print("calculating coefficients")
#alpha = ForwardSub(U,trainAEs)
#alpha = BackSub(np.transpose(U),alpha)

alpha = ForwardBackwardSub(np.transpose(U),trainAEs)

#print(len(alpha))

npred = len(prediction)
#npred = 10

L = np.zeros((len(training),npred))

print("building L")

for i in range(len(training)):
    for j in range(npred):
        L[i][j] = Kernel(training[i],prediction[j],Sigma,KernelType) #+ Lambda*KroD(i,j)
#        L[i][j] = Kernel(training[i],training[j],Sigma) #+ Lambda*KroD(i,j)

#estAEs = np.zeros(len(prediction))
estAEs = np.transpose(L).dot(alpha)

print(len(estAEs))

#print(estAEs,predictAEs[0:9])
#print(estAEs,trainAEs[0:10])

#error = MAD(estAEs,trainAEs[0:10])
error = MAD(estAEs,predictAEs)
print(error)


