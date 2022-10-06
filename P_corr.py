# -*- coding: utf-8 -*-
"""
Created on Wed May 18 15:13:17 2022

@author: berej
"""

import numpy as np
import os
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from mpmath import *

def P_mag(B):
    return B**2/(8*np.pi)

N = 49
P_gaz = []
z = []
write_list = []
B = []
z0 = 50e5
H = 100e5
mp.dps = 300

f=open('model-ProjetL3.txt','r')
text = f.readlines()
for i in range(1,N+1):
    line = text[i].split()
    P_gaz.insert(0,float(line[3]))
    z.insert(0,float(line[1]))
f.close()

f=open('P corrected.txt','w')

for p in range(len(P_gaz)):
    if z[p] == 100e5:
        P_ext = P_gaz[p]

for i in range(200,1300,100):
    P_corr = []
    for p in range(len(P_gaz)):
        pmag = P_mag(i)
        beta = (P_ext - pmag)/pmag 
        P_corr.append(P_gaz[p]/(1+1/beta))
# =============================================================================
#         if z[p] > z0:
#             P_corr.append(float(P_gaz[p]-P_mag(i*exp(-(z[p]-z0)/H))))
#         else:
#             P_corr.append(P_gaz[p]-P_mag(i))
# =============================================================================
    write_list.append(P_corr)
    B.append(i)
    
f.write('B(en G)\t\t')
for i in range(200,1300,100):
    if i == 1200:
        f.write(str(i)+'\n')
    else:
        f.write(str(i)+'\t\t')

f.write('z(en cm)\n')

for p in range(len(write_list[0])):
    f.write('{:.3e}\t'.format(z[p]))
    for i in range(len(write_list)):
        if i == len(write_list)-1:
            f.write('{:.3e}\n'.format(write_list[i][p]))
        else:
            f.write('{:.3e}\t'.format(write_list[i][p]))
    

f.close()

# =============================================================================
# [xx,yy] = np.meshgrid(z, B)
# 
# fig = plt.figure('P_mag as a function of B and z')
# ax = fig.gca(projection='3d')
# ax.plot_surface(xx, yy, np.array(write_list))
# =============================================================================


plt.figure('P_mag as a function of B and z')
plt.title('P_mag as a function of B and z')
for i in range(len(write_list)):
    plt.plot(z,write_list[i])
plt.show()

