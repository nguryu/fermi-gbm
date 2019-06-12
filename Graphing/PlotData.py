from numpy import *
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
# Lists to store whitened data.
white_real = []
white_res_real = []
white_im = []
white_res_im = []

# *** READ IN DATA *** #
with open('/Users/RonnyNguyen/Documents/Github/ldasoft/galactic_binaries/src/FisherGalaxy/Confusion_XAE_1.dat', 'r') as data_file:
    # data_slice = [next(data_file) for x in xrange(6000)]  # Read in file up to line N.
    # for line in data_file:
    for line in data_file.read().split("\n")[::100]:  # Read only every N lines.
        line = line.strip()  # Remove whitespace.
        if not line:  # Skip empty lines.
            continue
        if not line.startswith('#'):  # Skip comments.
            col = line.split(' ')  # Delimiter
            f_confusion.append(float(col[0]))
            x_noise.append(float(col[1]))
            x_conf.append(float(col[2]))
            ae_noise.append(float(col[3]))
            ae_conf.append(float(col[4]))

with open('/Users/RonnyNguyen/Documents/Github/ldasoft/galactic_binaries/src/FisherGalaxy/Galaxy_XAE.dat', 'r') as data_file2:
    for line in data_file2.read().split("\n")[::100]:
        line = line.strip()
        if not line:
            continue
        if not line.startswith('#'):
            col = line.split(' ')
            f_galaxy.append(float(col[0]))
            x_real.append(float(col[1]))
            x_im.append(float(col[2]))
            a_real.append(float(col[3]))
            a_im.append(float(col[4]))
            e_real.append(float(col[5]))
            e_im.append(float(col[6]))

with open('/Users/RonnyNguyen/Documents/Github/ldasoft/galactic_binaries/src/FisherGalaxy/Galaxy_XAE_R1.dat', 'r') as data_file3:
    for line in data_file3.read().split("\n")[::100]:
        line = line.strip()
        if not line:
            continue
        if not line.startswith('#'):
            col = line.split(' ')
            res_f_galaxy.append(float(col[0]))
            res_x_real.append(float(col[1]))
            res_x_im.append(float(col[2]))
            res_a_real.append(float(col[3]))
            res_a_im.append(float(col[4]))
            res_e_real.append(float(col[5]))
            res_e_im.append(float(col[6]))

# *** PERFORM CALCULATIONS TO PASS INTO FUNCTIONS *** #
# Take logarithm of data to plot.
for i in range(len(f_galaxy)):
    # Galaxy_XAE.dat
    if log10(f_galaxy[i]) >= -4:  # Assign cut-off frequency
        log_f_galaxy.append(log10(f_galaxy[i]))
        y_ax.append(log10(a_real[i]**2 + a_im[i]**2))  # Calculated y_ax axis
        a_real_cut.append(a_real[i])
        a_im_cut.append(a_im[i])
    # Galaxy_XAE_R1.dat
    if log10(res_f_galaxy[i]) >= -4:
        log_res_f_galaxy.append(log10(res_f_galaxy[i]))
        res_y.append(log10(res_a_real[i]**2 + res_a_im[i]**2))  # Calculated residual y_ax axis
        res_a_real_cut.append(res_a_real[i])
        res_a_im_cut.append(res_a_im[i])

# "Whiten" the real and imaginary amplitudes.
for i in range(len(f_confusion)):
    white_real.append(a_real_cut[i] / sqrt(ae_noise[i]))
    white_im.append(a_im_cut[i] / sqrt(ae_noise[i]))
    white_res_real.append(res_a_real_cut[i] / sqrt(ae_noise[i]))
    white_res_im.append(res_a_im_cut[i] / sqrt(ae_noise[i]))

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
    sub1 = plt.subplot(311)  # row i, column j, plot k
    #sub1.set_xlim([-4, -1.2])
    sub1.set_ylim([-50, -32])
    sub1.title.set_text('A vs. f')
    sub1.plot(log_f_galaxy, y_ax, color = 'gray', alpha = 1)

    # Subplot 2
    sub2 = plt.subplot(312)
    #sub2.set_xlim([-4, -1.2])
    sub2.set_ylim([-50, -32])
    sub2.title.set_text('Residual A vs. Residual f')
    sub2.set_ylabel(r'$\rm{log}(\rm{Re}(A)^2 + \rm{Im}(A)^2)$', visible = True)
    sub2.plot(log_res_f_galaxy, res_y, color = 'red', alpha = 0.1)

    # Subplot 3
    sub3 = plt.subplot(313)
    #sub3.set_xlim([-4, -1.2])
    sub3.set_ylim([-50, -32])
    sub3.title.set_text('Superimposed')
    sub3.set_xlabel(r'$\rm{log}(f)$', visible = True)
    sub3.plot(log_f_galaxy, y_ax, color = 'gray', alpha = 1)
    sub3.plot(log_res_f_galaxy, res_y, color = 'red', alpha = 0.1)

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
    sub1 = plt.subplot(311)  # row i, column j, plot k
    #sub1.set_xlim([-4, -1.2])
    #sub1.set_ylim([1.5, 2.5])
    sub1.title.set_text('A vs. f')
    sub1.plot(log_f_galaxy, white_y, color='gray', alpha=1)

    # Subplot 2
    sub2 = plt.subplot(312)
    #sub2.set_xlim([-4, -1.2])
    #sub2.set_ylim([1.5, 2.5])
    sub2.title.set_text('Residual A vs. Residual f')
    sub2.set_ylabel(r'$w(f)^2 = \rm{Re}(A)^2 + \rm{Im}(A)^2$', visible=True)
    #sub2.set_ylabel(r'$\rm{log}(w(f)^2) = \rm{log}(\rm{Re}(A)^2 + \rm{Im}(A)^2)$', visible=True)
    sub2.plot(log_res_f_galaxy, white_res_y, color='red', alpha=0.1)

    # Subplot 3
    sub3 = plt.subplot(313)
    #sub3.set_xlim([-4, -1.2])
    #sub3.set_ylim([1.5, 2.5])
    sub3.title.set_text('Superimposed')
    sub3.set_xlabel(r'$\rm{log}(f)$', visible=True)
    sub3.plot(log_f_galaxy, white_y, color='gray', alpha=1)
    sub3.plot(log_res_f_galaxy, white_res_y, color='red', alpha=0.1)

    # Output
    fig2.savefig('white_amp_vs_freq.eps', format='eps', bbox_inches='tight', pad_inches=0.02, dpi=1000)

def ifft_white(log_f_galaxy, log_res_f_galaxy, white_real, white_im, white_res_real, white_res_im):
    white_y = []
    white_res_y = []
    x_ax = []
    res_x = []
    f_galaxy = []
    res_f_galaxy = []

    # *** CREATE WHITENED Y-AXIS *** #
    for i in range(len(log_f_galaxy)):
        white_y.append(white_real[i]**2 + white_im[i]**2)
        white_res_y.append(white_res_real[i]**2 + white_res_im[i]**2)

    # *** INVERSE FOURIER TRANSFORM *** #
    # Return to regular frequency from the the logarithmic form.
    for i in range(len(log_f_galaxy)):
        f_galaxy.append(10 ** log_f_galaxy[i])
        res_f_galaxy.append(10 ** log_res_f_galaxy[i])
    # Inverse Fourier transform of frequency simply returns time.
    for i in range(len(f_galaxy)):
        x_ax.append(1/f_galaxy[i])
        res_x.append(1/res_f_galaxy[i])

    # Inverse Fourier transform the whitened amplitude and the frequency domain.
    ifft_white = fft.ifft(white_y)
    ifft_res_white = fft.ifft(white_res_y)
    # Take the logarithm for plotting visualization purposes.
    log_ifft_white = log10(fft.ifft(white_y))
    log_ifft_res_white = log10(fft.ifft(white_res_y))

    # *** FIGURE SET-UP *** #
    xr = 5.
    yr = 8.
    fig3 = plt.figure(figsize=(xr, yr))  # (width, height)

    # Subplot 1
    sub1 = plt.subplot(311)  # row i, column j, plot k
    #sub1.set_xlim([-4, -1.2])
    #sub1.set_ylim([-0.05, 0.05])
    sub1.title.set_text('A vs. t')
    sub1.plot(x_ax, ifft_white, color='gray', alpha=1)

    # Subplot 2
    sub2 = plt.subplot(312)
    #sub2.set_xlim([-4, -1.2])
    #sub2.set_ylim([-0.05, 0.05])
    sub2.title.set_text('Residual A vs. Residual t')
    sub2.set_ylabel(r'$w(t) = \mathcal{L}^{-1}(w(f))$', visible=True)
    #sub2.set_ylabel(r'$\rm{log}(w(t)) = \rm{log}(\mathcal{L}^{-1}(w(f)))$', visible=True)
    sub2.plot(res_x, ifft_res_white, color='red', alpha=0.1)

    # Subplot 3
    sub3 = plt.subplot(313)
    #sub3.set_xlim([-4, -1.2])
    #sub3.set_ylim([-0.05, 0.05])
    sub3.title.set_text('Superimposed')
    sub3.set_xlabel(r'$t$', visible=True)
    sub3.plot(x_ax, ifft_white, color='gray', alpha=1)
    sub3.plot(res_x, ifft_res_white, color='red', alpha=0.1)

    # Output
    fig3.savefig('ifft_white.eps', format='eps', bbox_inches='tight', pad_inches=0.02, dpi=1000)

# *** CALL FUNCTIONS *** #
# amp_vs_freq(log_f_galaxy, y_ax, log_res_f_galaxy, res_y)
white_amp_vs_freq(log_f_galaxy, log_res_f_galaxy, white_real, white_im, white_res_real, white_res_im)
#ifft_white(log_f_galaxy, log_res_f_galaxy, white_real, white_im, white_res_real, white_res_im)

end = time.time()
print "\nRuntime:", end - start, "s"