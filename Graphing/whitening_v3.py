#imports
import numpy as np
import matplotlib.pyplot as plt
from sys import getsizeof
import csv
import sys, os, math
import csv
import h5py

minf= 128. / 10000
maxf= 256. / 10000

#minf = 0.0001
#maxf = 0.0002

maxf=float(maxf)
minf=float(minf)

os.system("gcc -O2 -o bandpass_v2 ./bandpass_v2.c -lm")
os.system("./bandpass_v2 " + str(minf) + " " + str(maxf))

# Data
fig = plt.figure()
f2=open('./Galaxy_XAE_white.dat','r')
data=np.loadtxt(f2)
f=data[:,0]
re=data[:,1]
im=data[:,2]
df=(f[2]-f[1])
dt=1
Nt=np.int(np.round(1/(df*dt)))
fup=Nt/2*df
fBin=df
fSeries=np.arange(10**(minf), 10**(maxf), fBin)
white_0=np.zeros(len(fSeries))+1j*np.zeros(len(fSeries))
for i in range(0,len(fSeries),1):
    white_0[i]=complex(re[i],im[i])
white_td=np.fft.irfft(white_0,n=Nt)
yr_s=31457280
dt=1/(2*(Nt//2+1)*(f[2]-f[1]))
tup=Nt*dt/yr_s
timeBin=dt/yr_s
timeSeries=np.arange(0., tup, timeBin)
plt.plot(timeSeries[0::1000],white_td[0:len(timeSeries)][0::1000],'r')
plt.xlabel("Time (Yrs)",Fontsize=10)
plt.ylabel("Whitened data (red), Residual (black)",Fontsize=10)

# Residual
f2=open('./Galaxy_XAE_R1_white.dat','r')
data=np.loadtxt(f2)
f=data[:,0]
re=data[:,1]
im=data[:,2]
df=(f[2]-f[1])
dt=1
Nt=np.int(np.round(1/(df*dt)))
fup=Nt/2*df
fBin=df
fSeries=np.arange(10**(minf), 10**(maxf), fBin)
white_0_r=np.zeros(len(fSeries))+1j*np.zeros(len(fSeries))
for i in range(0,len(fSeries),1):
    white_0_r[i]=complex(re[i],im[i])
white_td_r=np.fft.irfft(white_0_r,n=Nt)
yr_s=31457280
dt=1/(2*(Nt//2+1)*(f[2]-f[1]))
tup=Nt*dt/yr_s
timeBin=dt/yr_s
timeSeries=np.arange(0., tup, timeBin)
plt.plot(timeSeries[0::1000],white_td_r[0:len(timeSeries)][0::1000],'k')
#plt.xlabel("Time (yr)",Fontsize=15)

fig.savefig('tukey8t.png', format='png', bbox_inches='tight', pad_inches=0.02, dpi=500)
