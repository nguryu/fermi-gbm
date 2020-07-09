import scipy
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from gbm.time import Met
from gbm.data.poshist import PosHist

# ================================================== #
# *** GBM/SGR CATALOGUE *** #
# ================================================== #
name = []
obsid = []
ra = []
dec = []
start_time = []
proc_time = []
xrt_exp = []
uvot_exp = []
bat_exp = []
archv_date = []

# Read in data.
with open('/Users/RonnyNguyen/Desktop/FermiGBMSummer2020/data_dir/BrowseTargets.txt', 'r') as data_file:
    for line in data_file:
        line = line.strip()  # Remove whitespace.
        if not line:  # Skip empty lines.
            continue
        if not line.startswith('#'):
            col = line.split('|')
            name.append(col[1])
            obsid.append(col[2])
            ra.append(col[3])
            dec.append(col[4])
            start_time.append(col[5])
            proc_time.append(col[6])
            xrt_exp.append(col[7])
            uvot_exp.append(col[8])
            bat_exp.append(col[9])
            archv_date.append(col[10])

# Search for SGR sources.
sgr_name = []
sgr_obsid = []
sgr_ra = []
sgr_dec = []

for i in range(len(name)):
    if name[i][0:3] == 'SGR':
        sgr_name.append(name[i])
        sgr_obsid.append(obsid[i])
        sgr_ra.append(ra[i])
        sgr_dec.append(dec[i])

# Convert RA and DEC to decimal.
sgr_ra_deg = []
sgr_dec_deg = []

for i in range(len(sgr_name)):
    ra_conv = (float(sgr_ra[i][0:2]) + float(sgr_ra[i][3:5])/60. + float(sgr_ra[i][6:11])/3600.)*360/24.
    dec_conv = float(sgr_dec[i][0:3]) + float(sgr_dec[i][4:6])/60. + float(sgr_dec[i][7:12])/3600.
    sgr_ra_deg.append(round(ra_conv, 2))
    sgr_dec_deg.append(round(dec_conv, 2))

# ================================================== #
# *** POSITION HISTORY FILE *** #
# ================================================== #
# Read in data.
poshist = PosHist.open("/Users/RonnyNguyen/Desktop/FermiGBMSummer2020/data_dir/glg_poshist_all_200619_v00.fit")

# Convert UTC time to mission elapsed time (MET, seconds) for Fermi.
time = "2020-06-19T15:20:36"
met = Met.from_iso(time).met

# Check if inside SAA +/- 10 sec pad.
# SAA = South Atlantic Anomaly where Fermi-GBM turns off due to the high particle flux.
in_saa = poshist.get_saa_passage(met) | poshist.get_saa_passage(met + 10) | poshist.get_saa_passage(met - 10)

# Check if true location is occulted.
occulted = poshist.location_visible(sgr_ra_deg, sgr_dec_deg, met) == False

if in_saa:
  print("Fermi GBM is in SAA +/- 10 sec")
elif occulted:
  print("Position not visible to Fermi GBM")
else:
  print("Source should be visible to Fermi GBM")
