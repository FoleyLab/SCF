#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr  4 14:25:38 2018

@author: jay
"""

import numpy as np
from numpy import linalg as LA
import math


def Map2EArray(length, dim, mat):
    Mmat = np.zeros((dim,dim,dim,dim))
    for m in range(0, length):
        i = int(mat[m][0]-1)
        j = int(mat[m][1]-1)
        k = int(mat[m][2]-1)
        l = int(mat[m][3]-1)
        val = mat[m][4]
        Mmat[i][j][k][l] = val
        Mmat[j][i][k][l] = val
        Mmat[i][j][l][k] = val
        Mmat[j][i][l][k] = val
        Mmat[k][l][i][j] = val
        Mmat[l][k][i][j] = val
        Mmat[k][l][j][i] = val
        Mmat[l][k][j][i] = val

    return Mmat

def UpdateF(dim, D, Hcore, EI):
    F = np.zeros_like(D)
    for m in range(0, dim):
        for n in range(0, dim):
            
            som = 0.
            
            for l in range(0, dim):
                for s in range(0, dim):
                    
                    som = som + D[l][s]*(2*EI[m][n][l][s] - EI[m][l][n][s])
                    
            
            F[m][n] = Hcore[m][n] + som
    
    return F
    
    
def Map1EArray(dim, mat):
    idx = 0
    Mmat = np.zeros((dim,dim))
    for i in range(0,dim):
       for j in range(i,dim):
          ival = int(mat[idx][0]-1)
          jval = int(mat[idx][1]-1)
          Mmat[ival][jval] = mat[idx][2]
          Mmat[jval][ival] = mat[idx][2]
          idx = idx + 1
    return Mmat

def BuildOrthoganlize(Smat):
  ## Diagonalize Total Hamiltonian, store vectors in v array
  SiSrt = np.zeros_like(Smat)
  Svals, Svecs = LA.eig(Smat)
  idx = Svals.argsort()[::1]
  Svals = Svals[idx]
  v = Svecs[:,idx]
  for i in range(0,len(Svals)):
    SiSrt[i][i] = Svals[i]**(-1/2.)

  S_int = np.dot(SiSrt,np.transpose(v))
  OM = np.dot(v, S_int)
  return OM

def BuildDensity(nocc, dim, C):
    D = np.zeros_like(C)
    for i in range(0,dim):
        for j in range(0, dim):
            som = 0
            for m in range(0, nocc):
                som = som + C[i][m]*C[j][m]
            D[i][j] = som
    return D

def RHF_Energy(dim, Hcore, D, F, Enuc):
    En = 0

    for i in range(0, dim):
        for j in range(0, dim):
            En = En + D[i][j]*(Hcore[i][j] + F[i][j])
    
    return (En + Enuc)
    
    
nocc = 5
Enuc = np.loadtxt("enuc.dat")
T    = np.loadtxt("T.dat")
V1   = np.loadtxt("V1.dat")
V2   = np.loadtxt("V2.dat")
S    = np.loadtxt("S.dat")

tei = int(len(V2))
dim = int(S[len(S)-1][0])

### store 1e data as arrays
Smat = Map1EArray(dim, S)
Tmat = Map1EArray(dim, T)
V1mat = Map1EArray(dim, V1)

### 2E array
TEls = Map2EArray(tei, dim, V2)
### define Hcore
Hcore = Tmat + V1mat

### Build S^{-1/2} matrix
OM = BuildOrthoganlize(Smat)

### Build Guess Fock Matrix as (S^{-1/2})^t Hcore S^{-1/2}
temp = np.dot(np.transpose(OM), Hcore)
F_init = np.dot(temp, OM)
#print(F_init)
### Diagonalize Fock matrix
eps0, C = LA.eig(F_init)
idx = eps0.argsort()[::1]
eps0 = eps0[idx]
Cp0 = C[:,idx]

#print(Cp0)
C0 = np.dot(OM, Cp0)

D0 = BuildDensity(nocc, dim, C0)
#print(D0)

E0 = RHF_Energy(dim, Hcore, D0, F_init, Enuc)

count=0
die=1
while(die):
    count = count + 1
    ### Update Fock Matrix
    Fp = UpdateF(dim, D0, Hcore, TEls)
    
    ### Transform fock matrix (S^{-1/2})^t Fp S^{-1/2}
    temp = np.dot(Fp,OM)
    F    = np.dot(np.transpose(OM), temp)

    ###Diagonalize F matrix - order in terms of increasing MO energies
    eps, Cp0 = LA.eig(F)
    idx = eps.argsort()[::1]
    eps = eps[idx]
    Cp = Cp0[:,idx]
    
    ### Get New MO coefficients
    C = np.dot(OM, Cp)
    ### Build Density
    Dnew = BuildDensity(nocc, dim, C)
    ### Compute New Energy
    Enew = RHF_Energy(dim, Hcore, Dnew, Fp, Enuc)

    ### Get Density Difference
    DD = Dnew-D0
    DDnorm = LA.norm(DD)
    
    ### Get energy difference
    deltaE = Enew-E0

    print("RHF Energy    Delta E     Delta D  ")
    print(Enew, deltaE, DDnorm)
    E0 = Enew
    D0 = Dnew
    
    if (count==34):
        die= 0

       
#def SortMatrix(mat):
#    for i in range(0,len(mat)):
        
