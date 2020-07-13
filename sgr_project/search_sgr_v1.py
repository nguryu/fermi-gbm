# This program takes the (RA, DEC) of an SGR recorded from the Swift catalogue and compares
# it to position history data from Fermi to see if the SGR was visible to Fermi on that day.

import os
from gbm.time import Met
from gbm.data.poshist import PosHist

# Inputs
path = '/Users/RonnyNguyen/Desktop/FermiGBMSummer2020/data_dir/'
catalogue_file = path + 'sgr_list.txt'

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
with open(catalogue_file, 'r') as data_file:
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

# Convert RA and DEC to decimal.
ra_deg = []
dec_deg = []
for i in range(len(name)):
    ra_conv = (float(ra[i][0:2]) + float(ra[i][3:5])/60. + float(ra[i][6:11])/3600.)*360/24.
    dec_conv = float(dec[i][0:3]) + float(dec[i][4:6])/60. + float(dec[i][7:12])/3600.
    ra_deg.append(round(ra_conv, 2))
    dec_deg.append(round(dec_conv, 2))

# ================================================== #
# *** FERMI POSITION HISTORY *** #
# ================================================== #
# Convert to correct time format for from_iso() function.
time_conv = []
for i in range(len(name)):
    time_conv.append(start_time[i][0:10] + 'T' + start_time[i][11:19])

# Check if source is visible to Fermi or not.
within_saa = []
not_visible = []
visible = []
for i in range(len(name)):
    if start_time[i][0:10] >= '2010-01-01':  # Position files before 2010 would give KeyErrors.
        yy = start_time[i][2:4]
        mm = start_time[i][5:7]
        dd = start_time[i][8:10]

        # Check if v00 and v01 file.
        if os.path.isfile(path+'glg_poshist_all_'+yy+mm+dd+'_v00.fit') == True:
            poshist_file = path+'glg_poshist_all_'+yy+mm+dd+'_v00.fit'
        else:
            poshist_file = path+'glg_poshist_all_'+yy+mm+dd+'_v01.fit'
            
        poshist = PosHist.open(poshist_file)

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
          within_saa.append(name[i].strip())  # strip() to remove extra white space.
        elif occulted:
          print("Position not visible to Fermi GBM")
          not_visible.append(name[i].strip())
        else:
          print("Source should be visible to Fermi GBM")
          visible.append(name[i].strip())

print("Within SAA:\t", within_saa)
print("Not visible:\t", not_visible)
print("Visible:\t", visible)
