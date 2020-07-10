import scipy
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from gbm.time import Met
from gbm.data.poshist import PosHist

search_date = '2020-06-19'

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
with open('/Users/RonnyNguyen/Desktop/FermiGBMSummer2020/data_dir/BrowseTargetsJune619.txt', 'r') as data_file:
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
sgr_start_time = []
sgr_ra = []
sgr_dec = []
for i in range(len(name)):
    if name[i][0:3] == 'SGR' and start_time[i][0:10] == search_date:
        sgr_name.append(name[i])
        sgr_obsid.append(obsid[i])
        sgr_start_time.append(start_time[i])
        sgr_ra.append(ra[i])
        sgr_dec.append(dec[i])

# Convert RA and DEC to decimal.
ra_deg = []
dec_deg = []
for i in range(len(sgr_name)):
    ra_conv = (float(sgr_ra[i][0:2]) + float(sgr_ra[i][3:5])/60. + float(sgr_ra[i][6:11])/3600.)*360/24.
    dec_conv = float(sgr_dec[i][0:3]) + float(sgr_dec[i][4:6])/60. + float(sgr_dec[i][7:12])/3600.
    ra_deg.append(round(ra_conv, 2))
    dec_deg.append(round(dec_conv, 2))

# ================================================== #
# *** POSITION HISTORY FILE *** #
# ================================================== #
# Read in data.
poshist = PosHist.open('/Users/RonnyNguyen/Desktop/FermiGBMSummer2020/data_dir/glg_poshist_all_200619_v00.fit')

# Convert to correct time format for from_iso() function.
time_conv = []
for i in range(len(sgr_start_time)):
    time_conv.append(sgr_start_time[i][0:10] + 'T' + sgr_start_time[i][11:19])

# Check if source is visible to Fermi or not.
within_saa = []
not_visible = []
visible = []
for i in range(len(sgr_name)):
    # Convert UTC time to mission elapsed time (MET in seconds) for Fermi.
    time = time_conv[i]
    met = Met.from_iso(time).met

    # Check if inside SAA +/- 10 sec pad.
    # SAA = South Atlantic Anomaly where Fermi-GBM turns off due to the high particle flux.
    in_saa = poshist.get_saa_passage(met) | poshist.get_saa_passage(met + 10) | poshist.get_saa_passage(met - 10)

    # Check if true location is occulted.
    occulted = poshist.location_visible(ra_deg[i], dec_deg[i], met) == False

    if in_saa:
      print("Fermi GBM is in SAA +/- 10 sec")
      within_saa.append(sgr_name[i].strip())  # strip() to remove extra white space.
    elif occulted:
      print("Position not visible to Fermi GBM")
      not_visible.append(sgr_name[i].strip())
    else:
      print("Source should be visible to Fermi GBM")
      visible.append(sgr_name[i].strip())

print("Within SAA:\t", within_saa)
print("Not visible:\t", not_visible)
print("Visible:\t", visible)
