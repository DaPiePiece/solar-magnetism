# -*- coding: utf-8 -*-
"""
Created on Thu May 19 11:31:21 2022

@author: berej
"""

from genericpath import isfile
from symtable import Symbol
from wsgiref import simple_server
import numpy as np
import os
import sympy
from matplotlib import pyplot as plt
from matplotlib import image
from astropy.io import fits
from mpmath import *
from natsort import natsorted

N=49
lamb = 630e-9
lamb_micro = 0.630
c = 3.00e8
h = 6.63e-34
nu = c/lamb
k = 1.38e-23
mp.dps = 300
lamb_1 = 1/lamb_micro
lamb_5 = 1/lamb_micro**5
mu = -0.998742

def trapint(y,x,a,b):
    intsum = 0
    for i in range(a,b):
        trapval = (x[i]-x[i-1])*0.5*(y[i]+y[i-1])
        intsum += trapval
    return intsum

def blamb(T, h=h,c=c,lamb=lamb,nu=nu, k=k):
    return ((2*h*nu**3)/c**2)*(1/(exp(h*nu/(k*T))-1))

wpath = os.path.join(os.getcwd(),'k')

files = natsorted([f for f in os.listdir(wpath) if os.path.isfile(os.path.join(wpath, f))])

z = []
k = [[] for i in range(len(files))]
taulist = [[] for i in range(len(files))]
B = []
I = []
T = []
G = []

for i in range(200,1300,100):
    G.append(i)

f=open('model-ProjetL3.txt','r')
text = f.readlines()
for i in range(1,N+1):
    line = text[i].split()
    z.insert(0,float(line[1]))
    T.insert(0,float(line[2]))
f.close()

wpath = os.path.join(os.getcwd(),'k')

files = natsorted([f for f in os.listdir(wpath) if os.path.isfile(os.path.join(wpath, f))])

for i in range(len(files)):
    f = open(os.path.join(wpath,files[i]),'r')
    text = f.readlines()
    for j in range(N):
        k[i].insert(0,float(text[j]))
    f.close()

plt.figure('opacity(B)')
plt.title('opacity(B)')
for i in range(len(k)):
    plt.plot(z,k[i])
plt.ylabel('k(z)')
plt.xlabel('z(km)')
plt.yscale('log')
plt.show()

for i in range(len(k)):
    for h in range(1,N+1):
        taulist[i].append(trapint(k[i],z,h,N))

plt.figure('tau(B)')
plt.title('tau(B)')
for i in range(len(taulist)):
    plt.plot(z,taulist[i])
plt.ylabel('tau(z)')
plt.xlabel('z(cm)')
plt.show()
        
for tau in taulist:
    intsum = 0
    for i in reversed(range(0,N-1)):
        trapval = 0.5*((tau[i]-tau[i+1])*(blamb(T[i])*exp(-1*tau[i]/mu)+blamb(T[i+1])*exp(-1*tau[i+1]/mu)))/mu
        intsum += trapval
    I.append(intsum)
    
plt.figure('I(B) theoretical')
plt.title('I(B) theoretical')
plt.plot(G,np.array(I)/max(I))
plt.xlabel('Parallel B (in Gauss)')
plt.ylabel('Normalised intensity')
plt.show()
