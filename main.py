from genericpath import isfile
from symtable import Symbol
from wsgiref import simple_server
import PIL
import numpy as np
import os
import sympy
from matplotlib import pyplot as plt
from matplotlib import image
from astropy.io import fits
from mpmath import *

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

a00 = 0.75267
a01 = -0.265577
a10 = 0.93874
a11 = 0.265577
a15 = -0.004095
a20 = -1.89287
a25 = 0.012582
a30 = 2.42234
a35 = -0.017117
a40 = -1.71150
a45 = 0.011977
a50 = 0.49062
a55 = -0.003347

A0 = a00+a01*lamb_1
A1 = a10+a11*lamb_1+a15*lamb_5
A2 = a20+a25*lamb_5
A3 = a30+a35*lamb_5
A4 = a40+a45*lamb_5
A5 = a50+a55*lamb_5
print(A0,A1,A2,A3,A4,A5)

def trapint(y,x,a,b):
    intsum = 0
    for i in range(a,b):
        trapval = (x[i]-x[i-1])*0.5*(y[i]+y[i-1])
        intsum += trapval
    return intsum

def blamb(T, h=h,c=c,lamb=lamb,nu=nu, k=k):
    return ((2*h*nu**3)/c**2)*(1/(exp(h*nu/(k*T))-1))

def poly(x):
    pol = []
    for i in range(len(x)):
        pol.append(A5*x[i]**5 + A4*x[i]**4 + A3*x[i]**3 + A2*x[i]**2 + A1*x[i] + A0)
    return pol

osn = os.name

path = os.getcwd()

z = []
k = []
tau = []
B = []
#neg_I = []
#pos_I = []
I = []
T = []
pos_noise = []
neg_noise = []

print('Data')
print('-Hinode 2007')
print('-Hinode 2014')
data_choice = input('Select data set: ')

print('Models')
print('-calm sun')
print('-sun spots')
model_choice = input('Select model: ')

if model_choice == 'calm sun':
    
    f=open(path+'/model-ProjetL3.txt','r')
    text = f.readlines()
    for i in range(1,N+1):
        line = text[i].split()
        z.insert(0,float(line[1]))
        k.insert(0,float(line[7]))
        T.insert(0,float(line[2]))
    f.close()
        
if model_choice == 'sun spots':
    
    f=open(path+'/model1007.dat','r')
    text=f.readlines()
    N = len(text)-2
    for i in range(0,N):
        line = text[i].split()
        z.insert(0,float(line[1])*10**5)
        T.insert(0,float(line[2]))
        k.insert(0,float(line[5]))
    
    f.close()
    
#print(z)
#print(k)
#print(T)

if data_choice == 'Hinode 2007':
    
    if osn == 'posix':
        path+='/data2007/'
    if osn == 'nt':
        path+='\\data2007\\'
else:
    
    if osn == 'posix':
        path+='/data2014/'
    if osn == 'nt':
        path+='\\data2014\\'
    
#print(path)
#print('Pillow version :',PIL.__version__)

files = [f for f in os.listdir(path) if os.path.isfile(os.path.join(path, f))]
#print(files)

xmax, ymax = (0,0)
dx, dy = (64,64)
angle = np.zeros(len(files))
radiance = np.zeros(len(files))
sigma = np.zeros(len(files))
pos_angle = []
pos_radiance = []
neg_angle = []
neg_radiance = []
pos_sigma = []
neg_sigma = []

for q in range(len(files)):

    imname = path+files[q]
    hdul = fits.open(imname)

    data = hdul[0].data
    xmax, ymax = data.shape
    # print(xmax, ymax)
    # print(xmax//2,ymax//2)

    #av = 0
    #for i in range(xmax//2-dx//2, xmax//2+dx//2):
    #    for j in range(ymax//2-dy//2,ymax//2+dy//2):
    #        av += data[i][j]
    #av = av/((dx*dy)+1)
    
    temp_sum = []
    av_sig = 0
    sig = 0
    dx = 63
    ystep = -1
    
    for ys in range(3):
        xc = 0
        for s in range(xmax//63):
            av = 0
            for i in range(xc, xc+dx):
                for j in range((ymax+dx*ystep)//2-dy//2,(ymax+dx*ystep)//2+dy//2):
                    av += data[i][j]
            av = av/((dx*dy)+1)
            temp_sum.append(av)
            xc = xc+dx
        ystep += 1
        
    for r in temp_sum:
        av_sig += r
    
    av_sig = av_sig / len(temp_sum)
    
    angle_value = float(files[q].lstrip('image-mu').rstrip('.fits').strip('='))
    radiance[q] = av_sig
    angle[q] = angle_value
    
    for r in temp_sum:
        sig += (r-av_sig)**2
        
    sigma[q] = sqrt(sig/(len(temp_sum)-1))
    
    hdul.close()

for i in range(len(angle)):
    if angle[i] < 0:
        neg_angle.append(-1*angle[i])
        neg_radiance.append(radiance[i])
        neg_sigma.append(sigma[i])
    else:
        pos_angle.append(angle[i])
        pos_radiance.append(radiance[i])
        pos_sigma.append(sigma[i])

for h in range(1,N+1):
    tau.append(trapint(k,z,h,N))

angle = neg_angle+pos_angle
angle.sort()

for mu in angle:
    intsum = 0
    for i in reversed(range(0,N-1)):
        trapval = 0.5*((tau[i]-tau[i+1])*(blamb(T[i])*exp(-1*tau[i]/mu)+blamb(T[i+1])*exp(-1*tau[i+1]/mu)))/mu
        intsum += trapval
    I.append(intsum)

# =============================================================================
# for mu in pos_angle:
#     intsum = 0
#     for i in reversed(range(0,N-1)):
#         trapval = 0.5*((tau[i]-tau[i+1])*(blamb(T[i])*exp(-1*tau[i]/mu)+blamb(T[i+1])*exp(-1*tau[i+1]/mu)))/mu
#         intsum += trapval
#     pos_I.append(intsum)
# =============================================================================

plt.figure()
plt.title('k(z) for lambda = 630nm')
plt.plot(z,k,c='b',label='z(k)')
plt.xlabel('Depth (cm)')
plt.ylabel('Absorption coefficient (cm-1)')
plt.yscale('log')
plt.show()
plt.figure()
plt.title('numerical solution of tau(z)')
plt.plot(z,tau,c='r',label='tau(z)')
plt.xlabel('Depth (cm)')
plt.ylabel('Optical depth')
plt.show()
plt.figure()
plt.title('T(z)')
plt.xlabel('Depth (cm)')
plt.ylabel('Temperature (K)')
plt.plot(z,T,c='r')
plt.show()

for i in range(len(pos_radiance)):
    pos_noise.append(float(sqrt(neg_radiance[i]*(dx*dy+1))/(neg_radiance[i]*(dx*dy+1))))

for i in range(len(neg_radiance)):
    neg_noise.append(float(sqrt(neg_radiance[i]*(dx*dy+1))/(neg_radiance[i]*(dx*dy+1))))
    
r_posmax = pos_radiance[len(pos_radiance)-1]
r_negmax = neg_radiance[len(neg_radiance)-1]
pos_radiance = pos_radiance/r_posmax
neg_radiance = neg_radiance/r_negmax

I_max = max(I)
print('sunspots model/calm sun model average intensity : '+'{:0.3e}'.format(float(I[-1])))

for i in range(len(I)):
    I[i]= I[i]/I_max
    
for i in range(len(pos_sigma)):
    pos_sigma[i] = pos_sigma[i]/r_posmax
    
for i in range(len(neg_sigma)):
    neg_sigma[i] = neg_sigma[i]/r_negmax
    
poly = poly(angle)

# =============================================================================
# print(pos_noise)
# print(neg_noise)
# 
# print(pos_sigma)
# print(neg_sigma)
# =============================================================================

plt.figure('Photon Noise')
plt.title('Photon noise according to mu angle '+data_choice)
plt.scatter(pos_angle,pos_noise,c='r',label='Positive mu angles')
plt.scatter(neg_angle,neg_noise,c='b',label='Negative mu anges')
plt.xlabel('mu cos(theta)')
plt.ylabel('sigma/<N>')
plt.legend()
plt.show()

plt.figure('Solar intensity according to mu angle')
plt.title('Averange intensity according to mu angles '+data_choice)
#plt.scatter(pos_angle, pos_radiance, c='r', label='fits data averages for positive angles')
plt.errorbar(pos_angle,pos_radiance,fmt='x',yerr=pos_sigma,c='r',label='fits data averages for positive angles')
#plt.plot(pos_angle,pos_I,c='r',label='I(tau,mu) for positive angles')
#plt.plot(pos_angle,pos_poly, c='purple', label='5th order polynomial fit for positive angles')
#plt.xlabel('mu (cos(theta))')
#plt.ylabel('Normalised radiance')
#plt.legend()
#plt.show()

#plt.figure('Solar radiance according to negative mu angle')
#plt.title('Averange radiance according to negative mu angles '+data_choice)
#plt.scatter(neg_angle, neg_radiance, c='b', label='fits data averages for negative angles')
plt.errorbar(neg_angle,neg_radiance,fmt='x',yerr=neg_sigma,c='b',label='fits data averages for negative angles')
plt.plot(angle,I,c='black',label='I(tau,mu)')
plt.plot(angle,poly, c='green', label='5th order polynomial fit')
plt.xlabel('mu (cos(theta))')
plt.ylabel('Normalised radiance')
plt.legend()
plt.show()