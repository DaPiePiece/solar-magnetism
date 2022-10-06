# -*- coding: utf-8 -*-
"""
Created on Tue May 10 10:26:59 2022

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
from mpl_toolkits.mplot3d import Axes3D

lamb = 630e-9
lamb_micro = 0.630
c = 3.00e8
h = 6.63e-34
nu = c/lamb
k = 1.38e-23
mp.dps = 300
lamb_1 = 1/lamb_micro
lamb_5 = 1/lamb_micro**5

def average(l):
    return sum(l)/len(l)

def path(decon,year,idstr=''):
    path = os.getcwd()
    if decon:
        path = os.path.join(path, 'decondata'+idstr+year)
    else:
        path = os.path.join(path, 'data'+idstr+year)
    return path

def trapint(y,x,a,b):
    intsum = 0
    for i in range(a,b):
        trapval = (x[i]-x[i-1])*0.5*(y[i]+y[i-1])
        intsum += trapval
    return intsum

def blamb(T, h=h,c=c,lamb=lamb,nu=nu, k=k):
    return ((2*h*nu**3)/c**2)*(1/(exp(h*nu/(k*T))-1))

def matrix_max(N):
    maxval = 0
    for i in N:
        for j in i:
            if j > maxval:
                maxval = j
    return maxval

intensity = []
bpl = []
iupfit = []
bupfit = []
idownfit = []
bdownfit = []
angle = []
avint = []
avinthighmag = []
P_corr = []
N = 49
z = []
k = [[] for i in range(11)]
taulist = [[] for i in range(11)]
B = []
I = []
T = []
G = [i for i in range(200,1300,100)]

print('Available data sets :')
print('2008')
print('2014')

year = input('select data set : ')

if year == '2014':
    alpha = 0.83 #calculating alpha dynamically is pointless in this application, we only have two sets of data.
else:
    alpha = 1

print('Deconvoluted?')

decon_choice = input('yes or no? (defaults to yes) : ')

if decon_choice == 'no':
    decon = False
else:
    decon = True   

wpath = path(decon,year)

files = [f for f in os.listdir(wpath) if os.path.isfile(os.path.join(wpath, f))]

for q in range(len(files)):
    if decon:    
        angle_value = float(files[q].lstrip('image-muD').rstrip('.fits').strip('='))
    else:
        angle_value = float(files[q].lstrip('image-mu').rstrip('.fits').strip('='))
    angle.append(angle_value)

print('possible angles : ', angle)
imagenum = input('select mu value (mu = cos(theta), closer to the center is closer to 1) : ')

mu = abs(float(imagenum))

if decon:
    current_file='imageD-mu='+imagenum+'.fits'
else:
    current_file='image-mu='+imagenum+'.fits'
    
readuplimit = input('B mode limit? (leave empty for 1250G) : ')

if readuplimit == '':
    uplimit = 1250
else:
    uplimit = int(readuplimit)
    
P_mag = uplimit**2/(8*np.pi)

write = input('write data? (defaults to no) : ')

if write == 'yes':
    writef = True
else:
    writef = False


imname = os.path.join(wpath, current_file)
hdul = fits.open(imname)
hdul.info()

idata = hdul[0].data
x_var, y_var = idata.shape

x_arr = [i for i in range(x_var)]
y_arr = [i for i in range(y_var)]

plt.figure('Heatmap of fits data for '+year+' '+current_file)
plt.title('Heatmap of fits data for '+year+' '+current_file)
plt.imshow(idata, cmap='hot', interpolation='nearest')
plt.colorbar()
plt.show()

hdul.close()

wpath = path(decon, year, 'mag')

if decon:
    current_file='BFe2D-mu='+imagenum+'.fits'
else:
    current_file='BFe2-mu='+imagenum+'.fits'

imname = os.path.join(wpath, current_file)
hdul = fits.open(imname)
hdul.info()

magdata = np.absolute(hdul[0].data)
xmax, ymax = magdata.shape
#print(xmax, ymax)
#print(xmax//2,ymax//2)

#print(magdata)
plt.figure('Heatmap of magnetic field data for '+year+' '+current_file)
plt.title('Heatmap of magnetic field data for '+year+' '+current_file)
plt.imshow(magdata, cmap='hot', interpolation='nearest')
plt.colorbar()
plt.show()

hdul.close()

for i in range(xmax):
    for j in range(ymax):
        
        mag = magdata[i][j]
        its = idata[i][j]*(1/alpha)
        
        if mag <= 100:
            avint.append(its)
        
        if mag >= 2000:
            avinthighmag.append(its)
        
        elif mag > 200 and mag <= uplimit:
            iupfit.append(its)
            bupfit.append(mag)
          
        elif mag > uplimit:
            idownfit.append(its)
            bdownfit.append(mag*(np.absolute(1/float(imagenum))))
        
        bpl.append(mag)
        intensity.append(its)

plt.figure('I(B) ' +current_file)
plt.title('Intensity according to parallel B for '+year+' '+current_file)
plt.scatter(bpl,intensity)
plt.xlabel('B parallel (in Gauss)')
plt.ylabel('I (in hits per pixel)')
plt.show()

# =============================================================================
# plt.figure('hist of B '+current_file)
# plt.title('Histogram of parallel B for '+year+' '+current_file)
# plt.yscale('log')
# plt.hist(bpl,bins=1000)
# plt.show()
# =============================================================================

# =============================================================================
# plt.figure('linear fit '+current_file)
# plt.title('Linear fit for intensity according to parallel B for '+year+' '+current_file)
# =============================================================================

(m1, b1), cov1 = np.polyfit(bupfit, iupfit, 1, cov=True)
u1 = np.sqrt(np.diag(cov1))[0]
x1 = np.sort(np.array(bupfit))
try:
    (m2, b2), cov2 = np.polyfit(bdownfit, idownfit, 1, cov=True)
    u2 = np.sqrt(np.diag(cov2))[0]
    x2 = np.sort(np.array(bdownfit))
    #plt.scatter(bdownfit, idownfit, color='lightsalmon', marker='o')
    #plt.plot(x2, m2*x2+b2, 'r', label = str(m2)+'*B+'+str(b2))
    rm2 = str(round(m2,3))
    ru2 = str(round(u2,3))
except:
    rm2 = 'N/A'
    ru2 = 'N/A'

y1 = m1*x1+b1
average_int = average(avint)
max_int = min(y1)

# =============================================================================
# for i in range(200,1300,100):
#     G.append(i)
# =============================================================================

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

for i in range(len(k)):
    for h in range(1,N+1):
        taulist[i].append(trapint(k[i],z,h,N))
        
for tau in taulist:
    intsum = 0
    for i in reversed(range(0,N-1)):
        trapval = 0.5*((tau[i]-tau[i+1])*(blamb(T[i])*exp(-1*tau[i]/mu)+blamb(T[i+1])*exp(-1*tau[i+1]/mu)))/mu
        intsum += trapval
    I.append(intsum)
    
plt.figure('I(B) theoretical')
plt.title('I(B) theoretical')
plt.scatter(bupfit, iupfit/max_int, color='lightblue', marker='o')
plt.plot(x1, y1/max_int, 'b', label = str(m1)+'*B+'+str(b1))
plt.plot(G,np.array(I)/min(I), label = 'theoretical I(B)')
plt.xlabel('Parallel B (in Gauss)')
plt.ylabel('Normalised intensity')
plt.legend()
plt.show()


#(os.getcwd())

if writef:
    f = open('results'+year+'.txt','a')
    f.write('mu = {mu}, upcoef = {m1:.3f}, u(upcoef) = {u1:.3f}, downcoef = {rm2}, u(downcoef) = {ru2} | avint = {avint:.3f}\n'.format(mu = imagenum, m1 = float(m1), u1 = float(u1), rm2 = rm2, ru2 = ru2, avint = float(average(avint))))
    f.close()
    
print('intensity average between -100 and 100 Gauss: '+str(average(avint)))
print('intesnity average higher than 2000 Gauss : '+str(average(avinthighmag)))
print('high/low field coefficient : '+str(average(avinthighmag)/average(avint)))

# =============================================================================
# [yy,xx] = np.meshgrid(y_arr, x_arr)
# 
# for i in range(len(idata)):
#     for j in range(len(idata[0])):
#         if idata[i][j] < 24000:
#             idata[i][j] = 24000
# 
# fig = plt.figure('magdata peaks and valleys')
# ax = fig.gca(projection='3d')
# ax.set_box_aspect([len(x_arr)/len(y_arr),1,1])  
# ax.plot_surface(xx, yy, idata/matrix_max(idata), cmap = 'hot')
# ax.plot_surface(xx, yy, magdata/matrix_max(magdata), cmap = 'gnuplot')
# =============================================================================
