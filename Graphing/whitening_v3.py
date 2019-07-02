import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from sys import getsizeof
import csv
import sys, os, math
import csv
import h5py
import time

start = time.time()

# Set min and max frequencies.
minf = 1. / 10000
maxf = 2. / 10000
outfile_name = "tukey_power1t.png"

# Run C programs to create tukey windows.
os.system("gcc -O2 -o bandpass_v2 ./bandpass_v2.c -lm")
os.system("./bandpass_v2 " + str(minf) + " " + str(maxf))

# *** DATA *** #
fields = ['f', 're', 'im']
df = pd.read_csv('./Galaxy_XAE_white.dat', skiprows=0, sep=" ", names=fields)
f = np.array(df.f)
re = np.array(df.re)
im = np.array(df.im)
# y-axis for power spectrum.
y_ax = np.abs(re**2 + im**2)
# Create frequency range.
df = (f[2]-f[1])
dt = 1
Nt = np.int(np.round(1/(df*dt)))
fup = Nt/2*df
fBin = df
fSeries = np.arange(10**(minf), 10**(maxf), fBin)
# Create the white y-axis as a complex array.
white_0 = np.zeros(len(fSeries))+1j*np.zeros(len(fSeries))
for i in range(0,len(fSeries), 1):
    white_0[i] = complex(re[i], im[i])
white_td = np.fft.irfft(white_0, n = Nt)
# Create time range.
yr_s = 31457280
dt = 1/(2*(Nt//2+1)*(f[2]-f[1]))
tup = Nt*dt/yr_s
timeBin = dt/yr_s
timeSeries = np.arange(0., tup, timeBin)
# Graphing.
xr = 8.
yr = 5.
fig = plt.figure(figsize = (xr, yr))  # (width, height)
#plt.plot(fSeries[0::1000], y_ax[0:len(fSeries)][0::1000], 'red')
#plt.ylabel("Whitened data (red), Residual (black)", fontsize = 10)
sub1 = plt.subplot(211)  # Row i, Column j, Plot k
sub1.set_ylabel(r'$\rm{Re}(A)^2 + \rm{Im}(A)^2$ Data', fontsize = 10)
sub1.plot(timeSeries[::1000], y_ax[0:len(timeSeries)][0::1000], color = "red")

# *** RESIDUAL *** #
fields2 = ['f2', 're2', 'im2']
df2 = pd.read_csv('./Galaxy_XAE_white.dat', skiprows=0, sep=" ", names=fields2)
f2 = np.array(df2.f2)
re2 = np.array(df2.re2)
im2 = np.array(df2.im2)
# y-axis for the power spectrum.
res_y = np.abs(re2**2 + im2**2)
# Get frequency range.
df = (f2[2]-f2[1])
dt = 1
Nt = np.int(np.round(1/(df*dt)))
fup = Nt/2*df
fBin = df
fSeries = np.arange(10**(minf), 10**(maxf), fBin)
# Create the white y-axis as a complex array.
white_0_r = np.zeros(len(fSeries))+1j*np.zeros(len(fSeries))
for i in range(0,len(fSeries), 1):
    white_0_r[i] = complex(re2[i], im2[i])
white_td_r=np.fft.irfft(white_0_r, n = Nt)
# Get time range.
yr_s = 31457280
dt = 1/(2*(Nt//2+1)*(f[2]-f[1]))
tup = Nt*dt/yr_s
timeBin = dt/yr_s
timeSeries = np.arange(0., tup, timeBin)
# Graphing.
#plt.plot(timeSeries[0::1000], white_td_r[0:len(timeSeries)][0::1000], 'k')
#plt.plot(fSeries[0::1000], res_y[0:len(fSeries)][0::1000], color = 'black', linestyle = ':', linewidth = 1)
sub2 = plt.subplot(212)
sub2.set_ylabel(r'$\rm{Re}(A)^2 + \rm{Im}(A)^2$ Residual', fontsize = 10)
sub2.set_xlabel("Time (Yrs)", fontsize = 10)
sub2.plot(timeSeries[::1000], res_y[0:len(timeSeries)][0::1000], color = "black")

end = time.time()
print("Runtime: " + "{:.2f}".format(end - start) + "s\n")

# *** OUTPUT *** #
fig.savefig(outfile_name, format='png', bbox_inches='tight', pad_inches=0.02, dpi=500)
