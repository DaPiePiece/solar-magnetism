import PIL
import numpy as np
import os
from matplotlib import pyplot as plt
from matplotlib import image
from astropy.io import fits

osn = os.name

path = os.getcwd()

if osn == 'posix':
    path+='/data2007/'
if osn == 'nt':
    path+='\\data2007\\'
    
print(path)
print('Pillow version :',PIL.__version__)

xmax, ymax = (0,0)
dx, dy = (2,2)
arrav = np.zeros(60)
print(arrav)
d_values = []

current_file='image-mu=0.75.fits'

imname = path+current_file
hdul = fits.open(imname)
#image = image.imread(path+'image-mu=0.5.fits')
# print(image.dtype)
# print(image.shape)

# plt.imshow(image)
# plt.show()
hdul.info()

data = hdul[0].data
xmax, ymax = data.shape
print(xmax, ymax)
print(xmax//2,ymax//2)

for q in range(len(arrav)):
    av = 0
    for i in range(xmax//2-dx//2, xmax//2+dx//2):
        for j in range(ymax//2-dy//2,ymax//2+dy//2):
            av += data[i][j]
    arrav[q] = av/((dx*dy)+1)
    d_values.append(dx)
    dx = dx+2
    dy = dy+2
#print(data[379][1023])
print(av)

print(data)
plt.figure('Heatmap of fits data for '+current_file)
plt.title('Heatmap of fits data for '+current_file)
plt.imshow(data, cmap='hot', interpolation='nearest')
#plt.matshow(data)
plt.show()
#plt.savefig('Heatmap of fits data for '+current_file+'.png', bbox_inches='tight')
plt.figure('sample size')
plt.plot(d_values,arrav)
plt.title('Averange intensity according to sample window size')
plt.xlabel('pixels')
plt.ylabel('number of hits')
plt.show()

hdul.close()