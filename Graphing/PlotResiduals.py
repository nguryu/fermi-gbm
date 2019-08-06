# ============================================================== #
# Ronny Nguyen
# Aug 6, 2019
# This code plots the whitened residual data for 2, 4, and 8 years
# of observation time.
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
fields1 = ['f_conf_2yr', 'x_noise_2yr', 'x_conf_2yr', 'ae_noise_2yr', 'ae_conf_2yr']
df1 = pd.read_csv('/Users/RonnyNguyen/Documents/Github/ldasoft/galactic_binaries/src/FisherGalaxy/2_year_data/Confusion_XAE_1.dat',
                 skiprows=0,
                 sep=" ",
                 names=fields1)

fields2 = ['res_f_galaxy_2yr', 'res_x_real_2yr', 'res_x_im_2yr', 'res_a_real_2yr', 'res_a_im_2yr', 'res_e_real_2yr', 'res_e_im_2yr']
df2 = pd.read_csv('/Users/RonnyNguyen/Documents/Github/ldasoft/galactic_binaries/src/FisherGalaxy/2_year_data/Galaxy_XAE_R1.dat',
                  skiprows=0,
                  sep=" ",
                  names=fields2)

fields3 = ['f_conf_4yr', 'x_noise_4yr', 'x_conf_4yr', 'ae_noise_4yr', 'ae_conf_4yr']
df3 = pd.read_csv('/Users/RonnyNguyen/Documents/Github/ldasoft/galactic_binaries/src/FisherGalaxy/4_year_data/Confusion_XAE_1.dat',
                 skiprows=0,
                 sep=" ",
                 names=fields3)

fields4 = ['res_f_galaxy_4yr', 'res_x_real_4yr', 'res_x_im_4yr', 'res_a_real_4yr', 'res_a_im_4yr', 'res_e_real_4yr', 'res_e_im_4yr']
df4 = pd.read_csv('/Users/RonnyNguyen/Documents/Github/ldasoft/galactic_binaries/src/FisherGalaxy/4_year_data/Galaxy_XAE_R1.dat',
                  skiprows=0,
                  sep=" ",
                  names=fields4)

fields5 = ['f_conf_8yr', 'x_noise_8yr', 'x_conf_8yr', 'ae_noise_8yr', 'ae_conf_8yr']
df5 = pd.read_csv('/Users/RonnyNguyen/Documents/Github/ldasoft/galactic_binaries/src/FisherGalaxy/8_year_data/Confusion_XAE_1.dat',
                 skiprows=0,
                 sep=" ",
                 names=fields5)

fields6 = ['res_f_galaxy_8yr', 'res_x_real_8yr', 'res_x_im_8yr', 'res_a_real_8yr', 'res_a_im_8yr', 'res_e_real_8yr', 'res_e_im_8yr']
df6 = pd.read_csv('/Users/RonnyNguyen/Documents/Github/ldasoft/galactic_binaries/src/FisherGalaxy/8_year_data/Galaxy_XAE_R1.dat',
                  skiprows=0,
                  sep=" ",
                  names=fields6)

# Put data into arrays.
f_conf_2yr = np.array(df1.f_conf_2yr)
x_noise_2yr = np.array(df1.x_noise_2yr)
x_conf_2yr = np.array(df1.x_conf_2yr)
ae_noise_2yr = np.array(df1.ae_noise_2yr)
ae_conf_2yr = np.array(df1.ae_conf_2yr)
res_f_galaxy_2yr = np.array(df2.res_f_galaxy_2yr)
res_x_real_2yr = np.array(df2.res_x_real_2yr)
res_x_im_2yr = np.array(df2.res_x_im_2yr)
res_a_real_2yr = np.array(df2.res_a_real_2yr)
res_a_im_2yr = np.array(df2.res_a_im_2yr)
res_e_real_2yr = np.array(df2.res_e_real_2yr)
res_e_im_2yr = np.array(df2.res_e_im_2yr)
f_conf_4yr = np.array(df3.f_conf_4yr)
x_noise_4yr = np.array(df3.x_noise_4yr)
x_conf_4yr = np.array(df3.x_conf_4yr)
ae_noise_4yr = np.array(df3.ae_noise_4yr)
ae_conf_4yr = np.array(df3.ae_conf_4yr)
res_f_galaxy_4yr = np.array(df4.res_f_galaxy_4yr)
res_x_real_4yr = np.array(df4.res_x_real_4yr)
res_x_im_4yr = np.array(df4.res_x_im_4yr)
res_a_real_4yr = np.array(df4.res_a_real_4yr)
res_a_im_4yr = np.array(df4.res_a_im_4yr)
res_e_real_4yr = np.array(df4.res_e_real_4yr)
res_e_im_4yr = np.array(df4.res_e_im_4yr)
f_conf_8yr = np.array(df5.f_conf_8yr)
x_noise_8yr = np.array(df5.x_noise_8yr)
x_conf_8yr = np.array(df5.x_conf_8yr)
ae_noise_8yr = np.array(df5.ae_noise_8yr)
ae_conf_8yr = np.array(df5.ae_conf_8yr)
res_f_galaxy_8yr = np.array(df6.res_f_galaxy_8yr)
res_x_real_8yr = np.array(df6.res_x_real_8yr)
res_x_im_8yr = np.array(df6.res_x_im_8yr)
res_a_real_8yr = np.array(df6.res_a_real_8yr)
res_a_im_8yr = np.array(df6.res_a_im_8yr)
res_e_real_8yr = np.array(df6.res_e_real_8yr)
res_e_im_8yr = np.array(df6.res_e_im_8yr)

end = time.time()
print("Checkpoint 1 " + "{:.2f}".format(end - start) + "s\n")

# Cast into real and imaginary components into a complex array.
a_res_2yr = np.zeros(len(res_f_galaxy_2yr)) + 1j*np.zeros(len(res_f_galaxy_2yr))
for i in range(len(res_f_galaxy_2yr)):
    a_res_2yr[i] = np.abs(complex(res_a_real_2yr[i]**2, res_a_im_2yr[i]**2))**2
    #a_res_2yr[i] = complex(res_a_real_2yr[i]**2, res_a_im_2yr[i]**2)

a_res_4yr = np.zeros(len(res_f_galaxy_4yr)) + 1j*np.zeros(len(res_f_galaxy_4yr))
for i in range(len(res_f_galaxy_4yr)):
    a_res_4yr[i] = np.abs(complex(res_a_real_4yr[i]**2, res_a_im_4yr[i]**2))**2
    #a_res_4yr[i] = complex(res_a_real_4yr[i]**2, res_a_im_4yr[i]**2)

a_res_8yr = np.zeros(len(res_f_galaxy_8yr)) + 1j*np.zeros(len(res_f_galaxy_8yr))
for i in range(len(res_f_galaxy_8yr)):
    a_res_8yr[i] = np.abs(complex(res_a_real_8yr[i]**2, res_a_im_8yr[i]**2))**2
    #a_res_8yr[i] = complex(res_a_real_8yr[i]**2, res_a_im_8yr[i]**2)


# ============================================================== #
# *** CALCULATIONS *** #
# ============================================================== #
lower = -2
upper = -4
f_conf1_2yr = []
a1_2yr = []
ae_noise1_2yr = []
ae_conf1_2yr = []
a_res1_2yr = []
f_conf1_4yr = []
a1_4yr = []
ae_noise1_4yr = []
ae_conf1_4yr = []
a_res1_4yr = []
f_conf1_8yr = []
a1_8yr = []
ae_noise1_8yr = []
ae_conf1_8yr = []
a_res1_8yr = []

# 2 Years
for i in range(len(f_conf_2yr)):
    if np.log10(f_conf_2yr[i]) <= lower and np.log10(f_conf_2yr[i]) >= upper:
        f_conf1_2yr.append(f_conf_2yr[i])
        ae_noise1_2yr.append(ae_noise_2yr[i])
        ae_conf1_2yr.append(ae_conf_2yr[i])
        a_res1_2yr.append(a_res_2yr[i])

# 4 Years
for i in range(len(f_conf_4yr)):
    if np.log10(f_conf_4yr[i]) <= lower and np.log10(f_conf_4yr[i]) >= upper:
        f_conf1_4yr.append(f_conf_4yr[i])
        ae_noise1_4yr.append(ae_noise_4yr[i])
        ae_conf1_4yr.append(ae_conf_4yr[i])
        a_res1_4yr.append(a_res_4yr[i])

# 8 Years
for i in range(len(f_conf_8yr)):
    if np.log10(f_conf_8yr[i]) <= lower and np.log10(f_conf_8yr[i]) >= upper:
        f_conf1_8yr.append(f_conf_8yr[i])
        ae_noise1_8yr.append(ae_noise_8yr[i])
        ae_conf1_8yr.append(ae_conf_8yr[i])
        a_res1_8yr.append(a_res_8yr[i])

# ============================================================== #
# *** WHITEN DATA *** #
# ============================================================== #
white_data_2yr = []
white_res_2yr = []
white_data_4yr = []
white_res_4yr = []
white_data_8yr = []
white_res_8yr = []

for i in range(len(f_conf1_2yr)):
    white_res_2yr.append(a_res_2yr[i] / np.sqrt(ae_noise1_2yr[i]))

for i in range(len(f_conf1_4yr)):
    white_res_4yr.append(a_res_4yr[i] / np.sqrt(ae_noise1_4yr[i]))

for i in range(len(f_conf1_8yr)):
    white_res_8yr.append(a_res_8yr[i] / np.sqrt(ae_noise1_8yr[i]))

end = time.time()
print("Checkpoint 2 " + "{:.2f}".format(end - start) + "s\n")

count_2yr = 0
count_4yr = 0
count_8yr = 0
cut_off = 0.0002

# Skip plotting first n elements in lists for better visualization.
# Comment this section out if you want to plot all points in array.
for i in range(len(f_conf1_2yr)):
    if f_conf1_2yr[i] <= cut_off:
        count_2yr += 1

for i in range(len(f_conf1_4yr)):
    if f_conf1_4yr[i] <= cut_off:
        count_4yr += 1

for i in range(len(f_conf1_8yr)):
    if f_conf1_8yr[i] <= cut_off:
        count_8yr += 1

# ============================================================== #
# *** PLOT RESIDUAL DATA *** #
# ============================================================== #
# Figure set-up.
xr = 8.
yr = 5.
fig1 = plt.figure(figsize=(xr, yr))  # (width, height)

# Plot data.
sub1 = plt.subplot(111)
sub1.title.set_text('Residual')
sub1.set_ylabel(r'$w(f)$', visible=True)
sub1.set_xlabel('Frequency (Hz)', visible=True)
sub1.loglog(f_conf1_2yr[count_2yr:][::100], white_res_2yr[count_2yr:][::100], color='red', zorder=0, label = '2 Years')
sub1.loglog(f_conf1_4yr[count_4yr:][::1000], white_res_4yr[count_4yr:][::1000], color='blue', zorder=1, label = '4 Years')
sub1.loglog(f_conf1_8yr[count_8yr:][::1000], white_res_8yr[count_8yr:][::1000], color='green', zorder=2, label = '8 Years')
pylab.legend(loc='upper right')

# Output
fig1.savefig('residual_log.png', format='png', bbox_inches='tight', pad_inches=0.02, dpi=500)

end = time.time()
print("Runtime: " + "{:.2f}".format(end - start) + "s\n")
