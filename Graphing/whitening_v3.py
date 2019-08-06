# ============================================================== #
# Ronny Nguyen
# Aug 6, 2019
# This code calls on bandpass_v2.c and creates a Tukey window for
# a desired min and max frequency range. It also does an inverse
# Fourier transform from frequency to time.
# ============================================================== #

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from sys import getsizeof
import sys, os, math
import csv
import h5py
from matplotlib.ticker import FormatStrFormatter
import time

start = time.time()

# Set min and max frequencies.
n = 5
minf = (2**1) / 10000
maxf = (2**n) / 10000
outfile_name = "tukey_2yr_" + str(n) + "t.png"  # Output file name.

# Run C programs to create Tukey windows.
os.system("gcc -O2 -o bandpass_v2 ./bandpass_v2.c -lm")
os.system("./bandpass_v2 " + str(minf) + " " + str(maxf))

# ============================================================== #
# *** DATA *** #
# ============================================================== #
fields = ['f', 're', 'im']
df = pd.read_csv('./Galaxy_XAE_white.dat', skiprows=0, sep=" ", names=fields)
f = np.array(df.f)
re = np.array(df.re)
im = np.array(df.im)

# Create frequency range.
df = f[2]-f[1]
dt = 1
Nt = np.int(np.round(1/(df*dt)))
fup = Nt/2*df
fBin = df
fSeries = np.arange(minf, maxf, fBin)

# Create the white y-axis as a complex array.
white_0 = np.zeros(len(fSeries))+1j*np.zeros(len(fSeries))
y_ax = np.zeros(len(fSeries))+np.zeros(len(fSeries))

for i in range(0, len(fSeries), 1):
    white_0[i] = np.abs(complex(re[i], im[i]))**2

# Take logarithm of axis.
log_white_0 = []
for i in range(0, len(white_0)):
    if white_0[i] > 0:
        log_white_0.append(np.log10(white_0[i]))
    else:
        log_white_0.append(0)

# Inverse Fourier transform.
ifft_white_y = np.fft.irfft(white_0, n = Nt)

# Create time range.
yr_s = 31457280. / 1
dt = 1/(2*(Nt/2+1)*(df))
tup = Nt*dt/yr_s
timeBin = dt/yr_s
timeSeries = np.arange(0., tup, timeBin)

# Graphing.
xr = 7.
yr = 5.
fig = plt.figure(figsize = (xr, yr))  # (width, height)
#plt.plot(tSeries[::1000], ifft_white_y[0:len(tSeries)][::1000], 'red')  # Plot inverse Fourier transform (time).
plt.plot(fSeries[::1000], log_white_0[0:len(fSeries)][::1000], 'red')  # Plot in terms of frequency.
plt.xlabel("Frequency (Hz)", fontsize = 10)
plt.ylabel("Whitened data (red), Residual (black)", fontsize = 10)

# ============================================================== #
# *** RESIDUAL *** #
# ============================================================== #
fields2 = ['f2', 're2', 'im2']
df2 = pd.read_csv('./Galaxy_XAE_R1_white.dat', skiprows=0, sep=" ", names=fields2)
f2 = np.array(df2.f2)
re2 = np.array(df2.re2)
im2 = np.array(df2.im2)

# Get frequency range.
df = f2[2] - f2[1]
dt = 1
Nt = np.int(np.round(1/(df*dt)))
fup = Nt/2*df
fBin = df
fSeries = np.arange(minf, maxf, fBin)

# Create the white y-axis as a complex array.
white_0_r = np.zeros(len(fSeries))+1j*np.zeros(len(fSeries))
res_y = np.zeros(len(fSeries))+np.zeros(len(fSeries))

for i in range(0,len(fSeries), 1):
    white_0_r[i] = np.abs(complex(re2[i], im2[i]))**2
    #white_0_r[i] = complex(re2[i], im2[i])

# Take logarithm of axis.
log_white_0_r = []
for i in range(0, len(white_0_r)):
    if white_0_r[i] > 0:
        log_white_0_r.append(np.log10(white_0_r[i]))
    else:
        log_white_0_r.append(0)

# Inverse Fourier transform.
ifft_white_res_y = np.fft.irfft(white_0_r, n = Nt)

# Get time range.
yr_s = 31457280. / 1
dt = 1/(2*(Nt/2+1)*(f[2]-f[1]))
tup = Nt*dt/yr_s
timeBin = dt / yr_s
timeSeries = np.arange(0., tup, timeBin)

# Graphing.
# plt.plot(tSeries[::1000], ifft_white_res_y[0:len(tSeries)][::1000], color = 'black')  # Plot inverse Fourier transform (time).
plt.plot(fSeries[::1000], log_white_0_r[0:len(fSeries)][::1000], color = 'black')  # Plot in terms of frequency.

# ============================================================== #
# *** OUTPUT *** #
# ============================================================== #
end = time.time()
print("Runtime: " + "{:.2f}".format(end - start) + "s\n")

fig.savefig(outfile_name, format='png', bbox_inches='tight', pad_inches=0.02, dpi=250)
