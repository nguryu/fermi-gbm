from numpy import *
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter, MultipleLocator
import time
start = time.time()

# *** INITIALIZE LISTS TO STORE DATA IN *** #
# Confusion_XAE_1.dat
f_confusion = []
x_noise = []
x_conf = []
ae_noise = []
ae_conf = []
# Galaxy_XAE.dat
f_galaxy = []
x_real = []
x_im = []
a_real = []
a_im = []
e_real = []
e_im = []
# Galaxy_XAE_R1.dat
res_f_galaxy = []
res_x_real = []
res_x_im = []
res_a_real = []
res_a_im = []
res_e_real = []
res_e_im = []
# Store calculations
y_ax = []
res_y = []
log_f_galaxy = []
log_res_f_galaxy = []
white_y = []
white_res_y = []
# These lists store the data after applying a cut-off frequency at 10**-4 Hz.
a_real_cut = []
res_a_real_cut = []
a_im_cut = []
res_a_im_cut = []
f_confusion_cut = []
ae_noise_cut = []
# Lists to store whitened data.
white_real = []
white_res_real = []
white_im = []
white_res_im = []

# *** READ IN DATA ** #
fields = ['f_confusion', 'x_noise', 'x_conf', 'ae_noise', 'ae_conf']
df = pd.read_csv('/Users/RonnyNguyen/Documents/Github/ldasoft/galactic_binaries/src/FisherGalaxy/Confusion_XAE_1.dat',
                 skiprows = 0,
                 sep = " ",
                 names = fields)

fields1 = ['f_galaxy', 'x_real', 'x_im', 'a_real', 'a_im', 'e_real', 'e_im']
df1 = pd.read_csv('/Users/RonnyNguyen/Documents/Github/ldasoft/galactic_binaries/src/FisherGalaxy/Galaxy_XAE.dat',
                 skiprows = 0,
                 sep = " ",
                 names = fields1)

fields2 = ['res_f_galaxy', 'res_x_real', 'res_x_im', 'res_a_real', 'res_a_im', 'res_e_real', 'res_e_im']
df2 = pd.read_csv('/Users/RonnyNguyen/Documents/Github/ldasoft/galactic_binaries/src/FisherGalaxy/Galaxy_XAE_R1.dat',
                 skiprows = 0,
                 sep = " ",
                 names = fields2)

# Put data into lists.
f_confusion = list(df.f_confusion)
x_noise = list(df.x_noise)
x_conf = list(df.x_conf)
ae_noise = list(df.ae_noise)
ae_conf = list(df.ae_conf)
f_galaxy = list(df1.f_galaxy)
x_real = list(df1.x_real)
x_im = list(df1.x_im)
a_real = list(df1.a_real)
a_im = list(df1.a_im)
e_real = list(df1.e_real)
e_im = list(df1.e_im)
res_f_galaxy = list(df2.res_f_galaxy)
res_x_real = list(df2.res_x_real)
res_x_im = list(df2.res_x_im)
res_a_real = list(df2.res_a_real)
res_a_im = list(df2.res_a_im)
res_e_real = list(df2.res_e_real)
res_e_im = list(df2.res_e_im)

# *** PERFORM CALCULATIONS TO PASS INTO FUNCTIONS *** #
# Count the number of frequencies cut off
num_zero = 0
res_num_zero = 0
conf_num_zero = 0

# Take logarithm of data to plot.
for i in range(len(f_galaxy)):
    # Galaxy_XAE.dat
    if log10(f_galaxy[i]) >= -4 and log10(f_galaxy[i]) <= -2:  # Assign cut-off frequency
        log_f_galaxy.append(log10(f_galaxy[i]))
        y_ax.append(log10(a_real[i]**2 + a_im[i]**2))  # Calculated y-axis
        a_real_cut.append(a_real[i])
        a_im_cut.append(a_im[i])
    else:
        num_zero += 1
    # Galaxy_XAE_R1.dat
    if log10(res_f_galaxy[i]) >= -4 and log10(res_f_galaxy[i]) <= -2:
        log_res_f_galaxy.append(log10(res_f_galaxy[i]))
        res_y.append(log10(res_a_real[i]**2 + res_a_im[i]**2))  # Calculated residual y-axis
        res_a_real_cut.append(res_a_real[i])
        res_a_im_cut.append(res_a_im[i])
    else:
        res_num_zero += 1

for i in range(len(f_confusion)-1):  # Ignore last last line in file.
    # Confusion_XAE.dat
    if log10(f_confusion[i]) >= -4 and log10(f_confusion[i]) <= -2:
        f_confusion_cut.append(f_confusion[i])
        ae_noise_cut.append(ae_noise[i])
    else:
        conf_num_zero += 1

print num_zero, res_num_zero, conf_num_zero
print len(a_real_cut), len(ae_noise_cut)

end = time.time()
print "\nRuntime:", end - start, "s"

# "Whiten" the real and imaginary amplitudes.
for i in range(len(log_f_galaxy)):
    # Galaxy_XAE.dat
    white_real.append(a_real_cut[i] / sqrt(ae_noise_cut[i]))
    white_im.append(a_im_cut[i] / sqrt(ae_noise_cut[i]))
    # Galaxy_XAE_R1.dat
    white_res_real.append(res_a_real_cut[i] / sqrt(ae_noise_cut[i]))
    white_res_im.append(res_a_im_cut[i] / sqrt(ae_noise_cut[i]))

# *** FIGURE PARAMETERS *** #
SMALL_SIZE = 8
MEDIUM_SIZE = 12
LARGE_SIZE = 16
plt.rc('font', size = SMALL_SIZE)          # Controls default text sizes
plt.rc('axes', titlesize = SMALL_SIZE)     # Font size of the axes title
plt.rc('axes', labelsize = SMALL_SIZE)     # Font size of the x and y labels
plt.rc('xtick', labelsize = SMALL_SIZE)    # Font size of the tick labels
plt.rc('ytick', labelsize = SMALL_SIZE)    # Font size of the tick labels
plt.rc('legend', fontsize = SMALL_SIZE)    # Legend font size
# Axes ticks
#majorFormatter = FormatStrFormatter('%1.1f')
#majorLocator = MultipleLocator(1000)  # Major tick intervals
#minorLocator = MultipleLocator(250)  # Minor tick intervals

def amp_vs_freq(log_f_galaxy, y_ax, log_res_f_galaxy, res_y):
    # *** FIGURE SET-UP *** #
    xr = 5.
    yr = 8.
    fig1 = plt.figure(figsize = (xr, yr)) # (width, height)

    # Subplot 1
    sub1 = plt.subplot(311)  # Row i, Column j, Plot K
    #sub1.set_xlim([-4, -1.2])
    sub1.set_ylim([-50, -32])
    sub1.title.set_text('A vs. f')
    sub1.plot(log_f_galaxy[::100], y_ax[::100], color = 'gray', alpha = 1)  # Plot only every N data points

    # Subplot 2
    sub2 = plt.subplot(312)
    #sub2.set_xlim([-4, -1.2])
    sub2.set_ylim([-50, -32])
    sub2.title.set_text('Residual A vs. Residual f')
    sub2.set_ylabel(r'$\rm{log}(\rm{Re}(A)^2 + \rm{Im}(A)^2)$', visible = True)
    sub2.plot(log_res_f_galaxy[::100], res_y[::100], color = 'red', alpha = 0.1)

    # Subplot 3
    sub3 = plt.subplot(313)
    #sub3.set_xlim([-4, -1.2])
    sub3.set_ylim([-50, -32])
    sub3.title.set_text('Superimposed')
    sub3.set_xlabel(r'$\rm{log}(f)$', visible = True)
    sub3.plot(log_f_galaxy[::100], y_ax[::100], color = 'gray', alpha = 1)
    sub3.plot(log_res_f_galaxy[::100], res_y[::100], color = 'red', alpha = 0.1)

    # Output
    fig1.savefig('amp_vs_freq.eps', format='eps', bbox_inches='tight', pad_inches=0.02, dpi=1000)

def white_amp_vs_freq(log_f_galaxy, log_res_f_galaxy, white_real, white_im, white_res_real, white_res_im):
    white_y = []
    white_res_y = []
    log_white_y = []
    log_white_res_y = []

    # *** CREATE Y-AXIS: Re(A)**2 + Im(A)**2 *** #
    for i in range(len(log_f_galaxy)):
        white_y.append(white_real[i]**2 + white_im[i]**2)
        white_res_y.append(white_res_real[i]**2 + white_res_im[i]**2)
    # Take the logarithm for plotting visualization purposes.
    for i in range(len(log_f_galaxy)):
        log_white_y.append(log10(white_y[i]))
        log_white_res_y.append(log10(white_res_y[i]))

    # *** FIGURE SET-UP *** #
    xr = 5.
    yr = 8.
    fig2 = plt.figure(figsize=(xr, yr))  # (width, height)
    # Subplot 1
    sub1 = plt.subplot(311)  # Row i, Column j, Plot K
    #sub1.set_xlim([-4, -1.2])
    #sub1.set_ylim([1.5, 2.5])
    sub1.title.set_text('A vs. f')
    sub1.plot(log_f_galaxy[::100], white_y[::100], color='gray', alpha=1)

    # Subplot 2
    sub2 = plt.subplot(312)
    #sub2.set_xlim([-4, -1.2])
    #sub2.set_ylim([1.5, 2.5])
    sub2.title.set_text('Residual A vs. Residual f')
    sub2.set_ylabel(r'$w(f)^2 = \rm{Re}(A)^2 + \rm{Im}(A)^2$', visible=True)
    #sub2.set_ylabel(r'$\rm{log}(w(f)^2) = \rm{log}(\rm{Re}(A)^2 + \rm{Im}(A)^2)$', visible=True)
    sub2.plot(log_res_f_galaxy[::100], white_res_y[::100], color='red', alpha=0.1)

    # Subplot 3
    sub3 = plt.subplot(313)
    #sub3.set_xlim([-4, -1.2])
    #sub3.set_ylim([1.5, 2.5])
    sub3.title.set_text('Superimposed')
    sub3.set_xlabel(r'$\rm{log}(f)$', visible=True)
    sub3.plot(log_f_galaxy[::100], white_y[::100], color='gray', alpha=1)
    sub3.plot(log_res_f_galaxy[::100], white_res_y[::100], color='red', alpha=0.1)

    # Output
    fig2.savefig('white_amp_vs_freq.eps', format='eps', bbox_inches='tight', pad_inches=0.02, dpi=1000)

def ifft_white(log_f_galaxy, log_res_f_galaxy, white_real, white_im, white_res_real, white_res_im,
               num_zero, res_num_zero):
    white_y = []
    white_res_y = []
    f_galaxy = []
    res_f_galaxy = []
    t_domain = []
    t_domain_res = []

    # *** CREATE WHITENED Y-AXIS *** #
    for i in range(len(log_f_galaxy)):
        white_y.append(complex(white_real[i], white_im[i]))
        white_res_y.append(complex(white_res_real[i], white_res_im[i]))

    # *** GET TIME DOMAIN FROM FREQUENCY *** #
    # Return to regular frequency from the the logarithmic form.
    for i in range(len(log_f_galaxy)):
        f_galaxy.append(10 ** log_f_galaxy[i])
        res_f_galaxy.append(10 ** log_res_f_galaxy[i])

    # Get time domain.
    df = f_galaxy[len(f_galaxy)-1] - f_galaxy[len(f_galaxy)-2]
    df_res = res_f_galaxy[len(res_f_galaxy)-1] - res_f_galaxy[len(res_f_galaxy)-2]
    dt = 1 / ((len(f_galaxy) + num_zero) * df)
    dt_res = 1 / ((len(res_f_galaxy) + res_num_zero) * df_res)

    for i in range(len(f_galaxy)):
        t_domain.append(i * dt)
        t_domain_res.append(i * dt_res)

    # *** INVERSE FOURIER TRANSFORM *** #
    ifft_white = fft.ifft(white_y)
    ifft_res_white = fft.ifft(white_res_y)

    # *** FIGURE SET-UP *** #
    xr = 5.
    yr = 8.
    fig3 = plt.figure(figsize=(xr, yr))  # (width, height)

    # Subplot 1
    sub1 = plt.subplot(311)  # Row i, Column j, Plot K
    sub1.title.set_text('A vs. t')
    sub1.plot(t_domain[::100], ifft_white[::100], color='gray', alpha=1)

    # Subplot 2
    sub2 = plt.subplot(312)
    sub2.title.set_text('Residual A vs. Residual t')
    sub2.set_ylabel(r'$w(t) = \mathcal{L}^{-1}(w(f))$', visible=True)
    #sub2.set_ylabel(r'$\rm{log}(w(t)) = \rm{log}(\mathcal{L}^{-1}(w(f)))$', visible=True)
    sub2.plot(t_domain_res[::100], ifft_res_white[::100], color='red', alpha=0.1)

    # Subplot 3
    sub3 = plt.subplot(313)
    sub3.title.set_text('Superimposed')
    sub3.set_xlabel(r'$t$', visible=True)
    sub3.plot(t_domain[::100], ifft_white[::100], color='gray', alpha=1)
    sub3.plot(t_domain_res[::100], ifft_res_white[::100], color='red', alpha=0.1)

    # Output
    fig3.savefig('ifft_white_cut.eps', format='eps', bbox_inches='tight', pad_inches=0.02, dpi=1000)

# *** CALL FUNCTIONS *** #
#amp_vs_freq(log_f_galaxy, y_ax, log_res_f_galaxy, res_y)
#white_amp_vs_freq(log_f_galaxy, log_res_f_galaxy, white_real, white_im, white_res_real, white_res_im)
ifft_white(log_f_galaxy, log_res_f_galaxy, white_real, white_im, white_res_real, white_res_im,
           num_zero, res_num_zero)

end = time.time()
print "\nRuntime:", end - start, "s"