# -*- coding: utf-8 -*-
"""
Created on Tue May 10 10:26:59 2022

@author: berej
"""

import PIL
import numpy as np
import os
from matplotlib import pyplot as plt
from matplotlib import image
from astropy.io import fits

def average(l):
    return sum(l)/len(l)

osn = os.name
xmax, ymax = (0,0)
dx, dy = (64,64)
intensity = []
bpl = []

iposupfit = []
bposupfit = []
iposdownfit = []
bposdownfit = []

inegupfit = []
bnegupfit = []
inegdownfit = []
bnegdownfit = []

angle = []
avint = []

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

path = os.getcwd()

if decon:
    if osn == 'posix':
        path+='/decondata'+year+'/'
    if osn == 'nt':
        path+='\\decondata'+year+'\\'
else:
    if osn == 'posix':
        path+='/data'+year+'/'
    if osn == 'nt':
        path+='\\data'+year+'\\'
    
files = [f for f in os.listdir(path) if os.path.isfile(os.path.join(path, f))]
for q in range(len(files)):
    if decon:    
        angle_value = float(files[q].lstrip('image-muD').rstrip('.fits').strip('='))
    else:
        angle_value = float(files[q].lstrip('image-mu').rstrip('.fits').strip('='))
    angle.append(angle_value)

print('possible angles : ', angle)
imagenum = input('select mu value (mu = cos(theta), closer to the center is closer to 1) : ')

if decon:
    current_file='imageD-mu='+imagenum+'.fits'
else:
    current_file='image-mu='+imagenum+'.fits'

imname = path+current_file
hdul = fits.open(imname)
hdul.info()

idata = hdul[0].data

plt.figure('Heatmap of fits data for '+year+' '+current_file)
plt.title('Heatmap of fits data for '+year+' '+current_file)
plt.imshow(idata, cmap='hot', interpolation='nearest')
plt.colorbar()
plt.show()

hdul.close()

path = os.getcwd()

if decon:
    if osn == 'posix':
        path+='/decondatamag'+year+'/'
    if osn == 'nt':
        path+='\\decondatamag'+year+'\\'
    current_file='BFe2D-mu='+imagenum+'.fits'
else:
    if osn == 'posix':
        path+='/datamag'+year+'/'
    if osn == 'nt':
        path+='\\datamag'+year+'\\'
    current_file='BFe2-mu='+imagenum+'.fits'

imname = path+current_file
hdul = fits.open(imname)
hdul.info()

magdata = hdul[0].data
xmax, ymax = magdata.shape
print(xmax, ymax)
print(xmax//2,ymax//2)

print(magdata)
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
        
        if mag > 300 and mag < 1400:
            iposupfit.append(its)
            bposupfit.append(mag)
            
        if mag < -300 and mag > -1400:
            inegupfit.append(its)
            bnegupfit.append(mag)
            
        if mag <= -1400:
            inegdownfit.append(its)
            bnegdownfit.append(mag)
         
        if mag >= 1400:
            iposdownfit.append(its)
            bposdownfit.append(mag)
        
        if mag >= -100 and mag <= 100:
            avint.append(its)
        
        bpl.append(mag)
        intensity.append(its)
        
(m1pos, b1pos), cov1pos = np.polyfit(bposupfit, iposupfit, 1, cov=True)
(m2pos, b2pos), cov2pos = np.polyfit(bposdownfit, iposdownfit, 1, cov=True)
(m1neg, b1neg), cov1neg = np.polyfit(bnegupfit, inegupfit, 1, cov=True)
if len(inegdownfit) > 10 :
    (m2neg, b2neg), cov2neg = np.polyfit(bnegdownfit, inegdownfit, 1, cov=True)

u1pos = np.sqrt(np.diag(cov1pos))
u2pos = np.sqrt(np.diag(cov2pos))
u1neg = np.sqrt(np.diag(cov1neg))
if len(inegdownfit) > 10 :
    u2neg = np.sqrt(np.diag(cov2neg))

x1pos = np.sort(np.array(bposupfit))
x2pos = np.sort(np.array(bposdownfit))
x1neg = np.sort(np.array(bnegupfit))
if len(inegdownfit) > 10 :
    x2neg = np.sort(np.array(bnegdownfit))
        
plt.figure('I(B) ' +current_file)
plt.title('Intensity according to parallel B for '+year+' '+current_file)
plt.scatter(bpl,intensity)
plt.xlabel('B parallel (in Gauss)')
plt.ylabel('I (in hits per pixel)')
plt.show()

plt.figure('linear fit '+current_file)
plt.title('Linear fit for intensity according to parallel B for '+year+' '+current_file)
plt.scatter(bposupfit, iposupfit, color='lightblue', marker='o')
plt.scatter(bnegupfit, inegupfit, color = 'lightblue', marker = 'o')
plt.plot(x1pos, m1pos*x1pos+b1pos, 'b', label = str(m1pos)+'*B+'+str(b1pos))
plt.plot(x1neg, m1neg*x1neg+b1neg, 'cyan', label = str(m1neg)+'*B+'+str(b1neg))
plt.scatter(bposdownfit, iposdownfit, color='lightsalmon', marker='o')
plt.plot(x2pos, m2pos*x2pos+b2pos, 'r', label = str(m2pos)+'*B+'+str(b2pos))
if len(inegdownfit) > 10 :
    plt.scatter(bnegdownfit, inegdownfit, color='lightsalmon', marker='o')
    plt.plot(x2neg, m2neg*x2neg+b2neg, 'magenta', label = str(m2neg)+'*B+'+str(b2neg))  

plt.xlabel('B parallel (in Gauss)')
plt.ylabel('I (in hits per pixel)')
plt.legend()
plt.show()

plt.figure('hist of B '+current_file)
plt.title('Histogram of parallel B for '+year+' '+current_file)
plt.yscale('log')
plt.hist(bpl,bins=1000)
plt.show()

print('intensity average between -100 and 100 Gauss: '+str(average(avint)))

if float(imagenum)<0:
    f = open('results south.txt','r')
    o_results = f.readlines()
    
    if len(inegdownfit) > 10 :
        if year == '2014':
            results = [o_results[0], o_results[1], '2014('+imagenum+'): upcoef = {m1pos:.3f}, u(upcoef) = {u1pos:0.3f}, downcoef = {m2pos:.3f}, u(downcoef) = {u2pos:.3f}\n'.format(m1pos=float(m1pos),u1pos=float(u1pos[0]),m2pos=float(m2pos),u2pos=float(u2pos[0])), o_results[3], o_results[4], '2014('+imagenum+'): upcoef = {m1neg:.3f}, u(upcoef) = {u1neg:0.3f}, downcoef = {m2neg:.3f}, u(downcoef) = {u2neg:.3f}\n'.format(m1neg=float(m1neg),u1neg=float(u1neg[0]),m2neg=float(m2neg),u2neg=float(u2neg[0]))]
        else:
            results = [o_results[0],'2008('+imagenum+'): upcoef = {m1pos:.3f}, u(upcoef) = {u1pos:0.3f}, downcoef = {m2pos:.3f}, u(downcoef) = {u2pos:.3f}\n'.format(m1pos=float(m1pos),u1pos=float(u1pos[0]),m2pos=float(m2pos),u2pos=float(u2pos[0])),o_results[2], o_results[3], '2008('+imagenum+'): upcoef = {m1neg:.3f}, u(upcoef) = {u1neg:0.3f}, downcoef = {m2neg:.3f}, u(downcoef) = {u2neg:.3f}\n'.format(m1neg=float(m1neg),u1neg=float(u1neg[0]),m2neg=float(m2neg),u2neg=float(u2neg[0])), o_results[5]]
    else:
        if year == '2014':
            results = [o_results[0], o_results[1], '2014('+imagenum+'): upcoef = {m1pos:.3f}, u(upcoef) = {u1pos:0.3f}, downcoef = {m2pos:.3f}, u(downcoef) = {u2pos:.3f}\n'.format(m1pos=float(m1pos),u1pos=float(u1pos[0]),m2pos=float(m2pos),u2pos=float(u2pos[0])), o_results[3], o_results[4], '2014('+imagenum+'): upcoef = {m1neg:.3f}, u(upcoef) = {u1neg:0.3f}, downcoef = {m2neg}, u(downcoef) = {u2neg}\n'.format(m1neg=float(m1neg),u1neg=float(u1neg[0]),m2neg='N/A',u2neg='N/A')]
        else:
            results = [o_results[0],'2008('+imagenum+'): upcoef = {m1pos:.3f}, u(upcoef) = {u1pos:0.3f}, downcoef = {m2pos:.3f}, u(downcoef) = {u2pos:.3f}\n'.format(m1pos=float(m1pos),u1pos=float(u1pos[0]),m2pos=float(m2pos),u2pos=float(u2pos[0])),o_results[2], o_results[3], '2008('+imagenum+'): upcoef = {m1neg:.3f}, u(upcoef) = {u1neg:0.3f}, downcoef = {m2neg}, u(downcoef) = {u2neg}\n'.format(m1neg=float(m1neg),u1neg=float(u1neg[0]),m2neg='N/A',u2neg='N/A'), o_results[5]]
        
    f.close()
    
    f=open('results south.txt','w')
    f.write(results[0]+results[1]+results[2]+results[3]+results[4]+results[5])
    f.close()
else:
    f = open('results north.txt','r')
    o_results = f.readlines()
    
    if len(inegdownfit) > 10 :
        if year == '2014':
            results = [o_results[0], o_results[1], '2014('+imagenum+'): upcoef = {m1pos:.3f}, u(upcoef) = {u1pos:0.3f}, downcoef = {m2pos:.3f}, u(downcoef) = {u2pos:.3f}\n'.format(m1pos=float(m1pos),u1pos=float(u1pos[0]),m2pos=float(m2pos),u2pos=float(u2pos[0])), o_results[3], o_results[4], '2014('+imagenum+'): upcoef = {m1neg:.3f}, u(upcoef) = {u1neg:0.3f}, downcoef = {m2neg:.3f}, u(downcoef) = {u2neg:.3f}\n'.format(m1neg=float(m1neg),u1neg=float(u1neg[0]),m2neg=float(m2neg),u2neg=float(u2neg[0]))]
        else:
            results = [o_results[0],'2008('+imagenum+'): upcoef = {m1pos:.3f}, u(upcoef) = {u1pos:0.3f}, downcoef = {m2pos:.3f}, u(downcoef) = {u2pos:.3f}\n'.format(m1pos=float(m1pos),u1pos=float(u1pos[0]),m2pos=float(m2pos),u2pos=float(u2pos[0])),o_results[2], o_results[3], '2008('+imagenum+'): upcoef = {m1neg:.3f}, u(upcoef) = {u1neg:0.3f}, downcoef = {m2neg:.3f}, u(downcoef) = {u2neg:.3f}\n'.format(m1neg=float(m1neg),u1neg=float(u1neg[0]),m2neg=float(m2neg),u2neg=float(u2neg[0])), o_results[5]]
    else:
        if year == '2014':
            results = [o_results[0], o_results[1], '2014('+imagenum+'): upcoef = {m1pos:.3f}, u(upcoef) = {u1pos:0.3f}, downcoef = {m2pos:.3f}, u(downcoef) = {u2pos:.3f}\n'.format(m1pos=float(m1pos),u1pos=float(u1pos[0]),m2pos=float(m2pos),u2pos=float(u2pos[0])), o_results[3], o_results[4], '2014('+imagenum+'): upcoef = {m1neg:.3f}, u(upcoef) = {u1neg:0.3f}, downcoef = {m2neg}, u(downcoef) = {u2neg}\n'.format(m1neg=float(m1neg),u1neg=float(u1neg[0]),m2neg='N/A',u2neg='N/A')]
        else:
            results = [o_results[0],'2008('+imagenum+'): upcoef = {m1pos:.3f}, u(upcoef) = {u1pos:0.3f}, downcoef = {m2pos:.3f}, u(downcoef) = {u2pos:.3f}\n'.format(m1pos=float(m1pos),u1pos=float(u1pos[0]),m2pos=float(m2pos),u2pos=float(u2pos[0])),o_results[2], o_results[3], '2008('+imagenum+'): upcoef = {m1neg:.3f}, u(upcoef) = {u1neg:0.3f}, downcoef = {m2neg}, u(downcoef) = {u2neg}\n'.format(m1neg=float(m1neg),u1neg=float(u1neg[0]),m2neg='N/A',u2neg='N/A'), o_results[5]]
        
    f.close()
    
    f=open('results north.txt','w')
    f.write(results[0]+results[1]+results[2]+results[3]+results[4]+results[5])
    f.close()