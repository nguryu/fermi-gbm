# ============================================================== #
# Ronny Nguyen
# Aug 6, 2019
# This code plots the whitened signal data as well as the residual
# data. Both are whitened in the process.
# ============================================================== #

import scipy
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter, MultipleLocator
from scipy import signal
import pylab
import time

start = time.time()

# ============================================================== #
# *** FIGURE PARAMETERS *** #
# ============================================================== #
SMALL_SIZE = 8
MEDIUM_SIZE = 12
LARGE_SIZE = 16
plt.rc('font', size=SMALL_SIZE)  # Controls default text sizes
plt.rc('axes', titlesize=SMALL_SIZE)  # Font size of the axes title
plt.rc('axes', labelsize=SMALL_SIZE)  # Font size of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE-2)  # Font size of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE-2)  # Font size of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE-2)  # Legend font size

# ============================================================== #
# *** READ IN DATA*** #
# ============================================================== #
fields = ['f_confusion', 'x_noise', 'x_conf', 'ae_noise', 'ae_conf']
df = pd.read_csv('/Users/RonnyNguyen/Documents/Github/ldasoft/galactic_binaries/src/FisherGalaxy/2_year_data/Confusion_XAE_1.dat',
                 skiprows=0,
                 sep=" ",
                 names=fields)

fields1 = ['f_galaxy', 'x_real', 'x_im', 'a_real', 'a_im', 'e_real', 'e_im']
df1 = pd.read_csv('/Users/RonnyNguyen/Documents/Github/ldasoft/galactic_binaries/src/FisherGalaxy/2_year_data/Galaxy_XAE.dat',
                  skiprows=0,
                  sep=" ",
                  names=fields1)

fields2 = ['res_f_galaxy', 'res_x_real', 'res_x_im', 'res_a_real', 'res_a_im', 'res_e_real', 'res_e_im']
df2 = pd.read_csv('/Users/RonnyNguyen/Documents/Github/ldasoft/galactic_binaries/src/FisherGalaxy/2_year_data/Galaxy_XAE_R1.dat',
                  skiprows=0,
                  sep=" ",
                  names=fields2)

# Put data into arrays.
f_confusion = np.array(df.f_confusion)
x_noise = np.array(df.x_noise)
x_conf = np.array(df.x_conf)
ae_noise = np.array(df.ae_noise)
ae_conf = np.array(df.ae_conf)
f_galaxy = np.array(df1.f_galaxy)
x_real = np.array(df1.x_real)
x_im = np.array(df1.x_im)
a_real = np.array(df1.a_real)
a_im = np.array(df1.a_im)
e_real = np.array(df1.e_real)
e_im = np.array(df1.e_im)
res_f_galaxy = np.array(df2.res_f_galaxy)
res_x_real = np.array(df2.res_x_real)
res_x_im = np.array(df2.res_x_im)
res_a_real = np.array(df2.res_a_real)
res_a_im = np.array(df2.res_a_im)
res_e_real = np.array(df2.res_e_real)
res_e_im = np.array(df2.res_e_im)

# Cast into real and imaginary components into a complex array.
a = np.zeros(len(f_galaxy)) + 1j*np.zeros(len(f_galaxy))
for i in range(len(f_galaxy)):
    a[i] = np.abs(complex(a_real[i], a_im[i]))**2

a_res = np.zeros(len(res_f_galaxy)) + 1j*np.zeros(len(res_f_galaxy))
for i in range(len(res_f_galaxy)):
    a_res[i] = np.abs(complex(res_a_real[i], res_a_im[i]))**2

# ============================================================== #
# *** CALCULATIONS *** #
# ============================================================== #
lower = -2
upper = -4
f_conf1 = []
a1 = []
ae_noise1 = []
ae_conf1 = []
a_res1 = []

for i in range(len(f_confusion)):
    if np.log10(f_confusion[i]) <= lower and np.log10(f_confusion[i]) >= upper:
        f_conf1.append(f_confusion[i])
        a1.append(a[i])
        ae_noise1.append(ae_noise[i])
        ae_conf1.append(ae_conf[i])
        a_res1.append(a_res[i])

# ============================================================== #
# *** WHITEN DATA *** #
# ============================================================== #
white_data = []
white_res = []

for i in range(len(f_conf1)):
    white_data.append(a1[i] / np.sqrt(ae_noise1[i]))
    white_res.append(a_res[i] / np.sqrt(ae_noise1[i]))

# ============================================================== #
# *** PLOT DATA *** #
# ============================================================== #
xr = 8.
yr = 5.
fig1 = plt.figure(figsize=(xr, yr))  # (width, height)

# Plot data.
sub1 = plt.subplot(111)
sub1.title.set_text('Whitened Data vs. Residual')
sub1.set_ylabel(r'$w(f)$', visible=True)
sub1.set_xlabel('Frequency (Hz)', visible=True)
sub1.loglog(f_conf1[1000:][::100], white_data[1000:][::100], color='red', zorder=0, label = 'Data')
sub1.loglog(f_conf1[1000:][::100], white_res[1000:][::100], color='black', zorder=1, label = 'Residual')
pylab.legend(loc='upper right')

# Output
fig1.savefig('white_2yr.png', format='png', bbox_inches='tight', pad_inches=0.02, dpi=500)

end = time.time()
print("Runtime: " + "{:.2f}".format(end - start) + "s\n")
