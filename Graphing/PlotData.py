import scipy
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter, MultipleLocator
from scipy import signal

# *** FIGURE PARAMETERS *** #
SMALL_SIZE = 8
MEDIUM_SIZE = 12
LARGE_SIZE = 16
plt.rc('font', size=SMALL_SIZE)  # Controls default text sizes
plt.rc('axes', titlesize=SMALL_SIZE)  # Font size of the axes title
plt.rc('axes', labelsize=SMALL_SIZE)  # Font size of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE-2)  # Font size of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE-2)  # Font size of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE-2)  # Legend font size

# Axes ticks
# majorFormatter = FormatStrFormatter('%1.1f')
# majorLocator = MultipleLocator(1000)  # Major tick intervals
# minorLocator = MultipleLocator(250)  # Minor tick intervals

# *** READ IN DATA ** #
fields = ['f_confusion', 'x_noise', 'x_conf', 'ae_noise', 'ae_conf']
df = pd.read_csv('/Users/RonnyNguyen/Documents/Github/ldasoft/galactic_binaries/src/FisherGalaxy/Confusion_XAE_1.dat',
                 skiprows=0,
                 sep=" ",
                 names=fields)

fields1 = ['f_galaxy', 'x_real', 'x_im', 'a_real', 'a_im', 'e_real', 'e_im']
df1 = pd.read_csv('/Users/RonnyNguyen/Documents/Github/ldasoft/galactic_binaries/src/FisherGalaxy/Galaxy_XAE.dat',
                  skiprows=0,
                  sep=" ",
                  names=fields1)

fields2 = ['res_f_galaxy', 'res_x_real', 'res_x_im', 'res_a_real', 'res_a_im', 'res_e_real', 'res_e_im']
df2 = pd.read_csv('/Users/RonnyNguyen/Documents/Github/ldasoft/galactic_binaries/src/FisherGalaxy/Galaxy_XAE_R1.dat',
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
    a[i] = complex(a_real[i], a_im[i])

a_res = np.zeros(len(res_f_galaxy)) + 1j*np.zeros(len(res_f_galaxy))
for i in range(len(res_f_galaxy)):
    a_res[i] = complex(res_a_real[i], res_a_im[i])

# *** PERFORM CALCULATIONS TO PASS INTO FUNCTIONS *** #
# Interpolation to align frequencies.
a_interp = np.interp(f_confusion, f_galaxy, a)
a_res_interp = np.interp(f_confusion, res_f_galaxy, a_res)
# Define desired frequency range.
f_conf1 = f_confusion[(np.log10(f_confusion) >= -2.2) & (np.log10(f_confusion) <= -1.9)]
a1 = a_interp[(np.log10(f_confusion) >= -2.2) & (np.log10(f_confusion) <= -1.9)]
ae_noise1 = ae_noise[(np.log10(f_confusion) >= -2.2) & (np.log10(f_confusion) <= -1.9)]
ae_conf1 = ae_conf[(np.log10(f_confusion) >= -2.2) & (np.log10(f_confusion) <= -1.9)]
a_res1 = a_res_interp[(np.log10(f_confusion) >= -2.2) & (np.log10(f_confusion) <= -1.9)]

print(len(a_res1))
print(len(f_conf1))

# *** WHITEN OUT DATA *** #
white = a1 / np.sqrt(ae_noise1)
white_res = a_res1 / np.sqrt(ae_noise1)
# Zero pad the lists for FFT.
z1 = len(f_confusion) - len(f_confusion[np.log10(f_confusion) <= -1.9])  # Get np.zeros after
z2 = len(f_confusion) - len(f_confusion[np.log10(f_confusion) >= -2.2])  # Get np.zeros before
white = np.hstack([np.zeros(z2), white])
white = np.hstack([white, np.zeros(z1 - 1)])
white_res = np.hstack([np.zeros(z2), white_res])
white_res = np.hstack([white_res, np.zeros(z1 - 1)])

def amp_vs_freq(log_f_galaxy, y_ax, log_res_f_galaxy, res_y):
    # *** FIGURE SET-UP *** #
    xr = 5.
    yr = 8.
    fig1 = plt.figure(figsize=(xr, yr))  # (width, height)

    # Subplot 1
    sub1 = plt.subplot(311)  # Row i, Column j, Plot K
    # sub1.set_xlim([-4, -1.2])
    sub1.set_ylim([-50, -32])
    sub1.title.set_text('A vs. f')
    sub1.plot(log_f_galaxy[::1], y_ax[::1], color='gray', alpha=1)  # Plot only every N data points

    # Subplot 2
    sub2 = plt.subplot(312)
    # sub2.set_xlim([-4, -1.2])
    sub2.set_ylim([-50, -32])
    sub2.title.set_text('Residual A vs. Residual f')
    sub2.set_ylabel(r'$\rm{log}(\rm{Re}(A)^2 + \rm{Im}(A)^2)$', visible=True)
    sub2.plot(log_res_f_galaxy[::1], res_y[::1], color='red', alpha=0.1)

    # Subplot 3
    sub3 = plt.subplot(313)
    # sub3.set_xlim([-4, -1.2])
    sub3.set_ylim([-50, -32])
    sub3.title.set_text('Superimposed')
    sub3.set_xlabel(r'$\rm{log}(f)$', visible=True)
    sub3.plot(log_f_galaxy[::1], y_ax[::1], color='gray', alpha=1)
    sub3.plot(log_res_f_galaxy[::1], res_y[::1], color='red', alpha=0.1)

    # Output
    fig1.savefig('amp_vs_freq.eps', format='eps', bbox_inches='tight', pad_inches=0.02, dpi=1000)

def white_amp_vs_freq(f_confusion, f_conf1, f_galaxy, res_f_galaxy, white, white_res):
    # *** FREQUENCY DOMAIN *** #
    df = f_confusion[2] - f_confusion[1]
    dt = 1  # Sampling time
    nt = np.int(np.round(1 / (df * dt)))
    f_up = nt / 2 * df
    f_bin = df
    f_series = np.arange(0., f_up, f_bin)

    # *** WHITE SPECTRA ** #
    tmp=np.zeros(2**(int(np.log2(nt)))-len(white))+1j*np.zeros(2**(int(np.log2(nt)))-len(white))
    white_y = np.hstack([white, tmp])
    tmp=np.zeros(2**(int(np.log2(nt)))-len(white_res))+1j*np.zeros(2**(int(np.log2(nt)))-len(white_res))
    white_res_y = np.hstack([white_res, tmp])

    f_series = f_series * 1e3
    white_y = np.absolute(white_y[0:len(f_series)])**2
    white_res_y = np.absolute(white_res_y[0:len(f_series)])**2

    # *** FIGURE SET-UP *** #
    xr = 5.
    yr = 8.
    fig2 = plt.figure(figsize=(xr, yr))  # (width, height)

    # Subplot 1
    sub1 = plt.subplot(311)  # Row i, Column j, Plot k
    sub1.title.set_text('A vs. f')
    sub1.loglog(f_series[::100], white_y[::100], color='gray')

    # Subplot 2
    sub2 = plt.subplot(312)
    sub2.title.set_text('Residual A vs. Residual t')
    sub2.set_ylabel(r'$w(f)^2 = \rm{Re}(A)^2 + \rm{Im}(A)^2$', visible=True)
    sub2.loglog(f_series[::100], white_res_y[::100], color='red')

    # Subplot 3
    sub3 = plt.subplot(313)
    sub3.title.set_text('Superimposed')
    sub3.set_xlabel(r'$f$', visible=True)
    sub3.loglog(f_series[::100], white_y[::100], color='gray')
    sub3.loglog(f_series[::100], white_res_y[::100], color='red')

    # Output
    fig2.savefig('white.eps', format='eps', bbox_inches='tight', pad_inches=0.02, dpi=500)

def ifft_white(f_confusion, f_conf1, f_galaxy, res_f_galaxy, white, white_res):
    # *** TIME AND FREQUENCY DOMAIN *** #
    # Get the frequency domain.
    df = f_confusion[2] - f_confusion[1]
    dt = 1  # Sampling time
    nt = np.int(np.round(1 / (df * dt)))
    f_up = nt / 2 * df
    f_bin = df
    f_series = np.arange(0., f_up, f_bin)
    # Get the time domain.
    yr_s = 31457280
    dt = 1 / (2*(nt//2+1) * (f_conf1[2] - f_conf1[1]))
    t_up = nt * dt / yr_s
    t_bin = dt / yr_s
    t_series = np.arange(0., t_up, t_bin)

    # *** WHITE SPECTRA ** #
    tmp=np.zeros(2**(int(np.log2(nt)))-len(white))+1j*np.zeros(2**(int(np.log2(nt)))-len(white))
    white_y = np.hstack([white, tmp])
    tmp=np.zeros(2**(int(np.log2(nt)))-len(white_res))+1j*np.zeros(2**(int(np.log2(nt)))-len(white_res))
    white_res_y = np.hstack([white_res, tmp])

    # *** INVERSE FOURIER TRANSFORM *** #
    ifft_white = np.fft.irfft(white_y, n = nt)
    ifft_res_white = np.fft.irfft(white_res_y, n = nt)
    ifft_white = ifft_white[0:len(t_series)]
    ifft_res_white = ifft_res_white[0:len(t_series)]

    # *** FIGURE SET-UP *** #
    xr = 5.
    yr = 8.
    fig3 = plt.figure(figsize=(xr, yr))  # (width, height)

    # Subplot 1
    sub1 = plt.subplot(311)  # Row i, Column j, Plot k
    sub1.title.set_text('A vs. t')
    sub1.plot(t_series[::100], ifft_white[::100], color='gray')

    # Subplot 2
    sub2 = plt.subplot(312)
    sub2.title.set_text('Residual A vs. Residual t')
    sub2.set_ylabel(r'$w(t) = \mathcal{L}^{-1}(w(f))$', visible=True)
    sub2.plot(t_series[::100], ifft_res_white[::100], color='red')

    # Subplot 3
    sub3 = plt.subplot(313)
    sub3.title.set_text('Superimposed')
    sub3.set_xlabel(r'$t$ (years)', visible=True)
    sub3.plot(t_series[::100], ifft_white[::100], color='gray')
    sub3.plot(t_series[::100], ifft_res_white[::100], color='red')

    # Output
    fig3.savefig('ifft_white_2mHz.eps', format='eps', bbox_inches='tight', pad_inches=0.02, dpi=500)

def tukey(f_confusion, f_conf1, f_galaxy, res_f_galaxy, white, white_res):
    # *** TIME AND FREQUENCY DOMAIN *** #
    # Get the frequency domain.
    df = f_confusion[2] - f_confusion[1]
    dt = 1  # Sampling time
    nt = np.int(np.round(1 / (df * dt)))
    f_up = nt / 2 * df
    f_bin = df
    f_series = np.arange(0., f_up, f_bin)

    # *** WHITE SPECTRA *** #
    '''
    tmp=np.zeros(2**(int(np.log2(nt)))-len(white))+1j*np.zeros(2**(int(np.log2(nt)))-len(white))
    white_y = np.hstack([white, tmp])
    tmp=np.zeros(2**(int(np.log2(nt)))-len(white_res))+1j*np.zeros(2**(int(np.log2(nt)))-len(white_res))
    white_res_y = np.hstack([white_res, tmp])
    '''

    # Remove zeros from data.
    white_y = [white[i] for i in range(len(white)) if white[i] != 0]
    white_res_y = [white_res[i] for i in range(len(white_res)) if white_res[i] != 0]

    print(len(white_y), len(white_res_y), len(f_conf1))

    # *** INVERSE FOURIER TRANSFORM *** #
    ifft_white = np.fft.irfft(white_y, n = nt)
    ifft_res_white = np.fft.irfft(white_res_y, n = nt)
    ifft_white = np.fft.fftshift(ifft_white)
    ifft_res_white = np.fft.fftshift(ifft_res_white)
    #fft.fftfreq

    # *** CREATE TUKEY WINDOW *** #
    window = signal.tukey(len(f_conf1))

    # *** FIGURE SET-UP *** #
    xr = 5.
    yr = 8.
    fig4 = plt.figure()  # (width, height)

    # Subplot 1
    sub1 = plt.subplot(111)  # Row i, Column j, Plot k
    sub1.title.set_text('Tukey')
    #sub1.set_ylim([-0.00001, 0.00001])
    #sub1.plot(window)
    sub1.plot(f_conf1, ifft_white[0:len(f_conf1)], color='gray')

    # Output
    fig4.savefig('tukey_window7.eps', format='eps', bbox_inches='tight', pad_inches=0.02, dpi=500)

# *** CALL FUNCTIONS *** #
#amp_vs_freq(log_f_galaxy, y_ax, log_res_f_galaxy, res_y)
#white_amp_vs_freq(f_confusion, f_conf1, f_galaxy, res_f_galaxy, white, white_res)
#ifft_white(f_confusion, f_conf1, f_galaxy, res_f_galaxy, white, white_res)
tukey(f_confusion, f_conf1, f_galaxy, res_f_galaxy, white, white_res)
