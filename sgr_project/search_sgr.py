# This program takes the (RA, DEC) of an SGR recorded from the Swift catalogue and compares
# it to position history data from Fermi to see if the SGR was visible to Fermi on that day.

import os
from gbm.time import Met
from gbm.data.poshist import PosHist

# Inputs
path = '/Users/RonnyNguyen/Desktop/FermiGBMSummer2020/data_dir/'
catalogue_file = path + 'swift_sgr_list.txt'

# ================================================== #
# *** SGR LIST *** #
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
within_saa_name = []
within_saa_obsid = []
within_saa_ra = []
within_saa_dec = []
within_saa_start_time = []
within_saa_proc_time = []
within_saa_xrt_exp = []
within_saa_uvot_exp = []
within_saa_bat_exp = []
within_saa_archv_date = []
not_visible_name = []
not_visible_obsid = []
not_visible_ra = []
not_visible_dec = []
not_visible_start_time = []
not_visible_proc_time = []
not_visible_xrt_exp = []
not_visible_uvot_exp = []
not_visible_bat_exp = []
not_visible_archv_date = []
visible_name = []
visible_obsid = []
visible_ra = []
visible_dec = []
visible_start_time = []
visible_proc_time = []
visible_xrt_exp = []
visible_uvot_exp = []
visible_bat_exp = []
visible_archv_date = []

for i in range(len(name)):
    if start_time[i][0:10] >= '2009-06-01':  # Position files before mid-2009 would give KeyErrors.
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
          # print("Fermi GBM is in SAA +/- 10 sec")
          within_saa_name.append(name[i])
          within_saa_obsid.append(obsid[i])
          within_saa_ra.append(ra[i])
          within_saa_dec.append(dec[i])
          within_saa_start_time.append(start_time[i])
          within_saa_proc_time.append(proc_time[i])
          within_saa_xrt_exp.append(xrt_exp[i])
          within_saa_uvot_exp.append(uvot_exp[i])
          within_saa_bat_exp.append(bat_exp[i])
          within_saa_archv_date.append(archv_date[i])
        elif occulted:
          # print("Position not visible to Fermi GBM")
          not_visible_name.append(name[i])
          not_visible_obsid.append(obsid[i])
          not_visible_ra.append(ra[i])
          not_visible_dec.append(dec[i])
          not_visible_start_time.append(start_time[i])
          not_visible_proc_time.append(proc_time[i])
          not_visible_xrt_exp.append(xrt_exp[i])
          not_visible_uvot_exp.append(uvot_exp[i])
          not_visible_bat_exp.append(bat_exp[i])
          not_visible_archv_date.append(archv_date[i])
        else:
          # print("Source should be visible to Fermi GBM")
          visible_name.append(name[i])
          visible_obsid.append(obsid[i])
          visible_ra.append(ra[i])
          visible_dec.append(dec[i])
          visible_start_time.append(start_time[i])
          visible_proc_time.append(proc_time[i])
          visible_xrt_exp.append(xrt_exp[i])
          visible_uvot_exp.append(uvot_exp[i])
          visible_bat_exp.append(bat_exp[i])
          visible_archv_date.append(archv_date[i])

# Write relevant data into new file
with open(path+'swift_within_saa.txt', 'w') as outfile1:
    datalist1 = [within_saa_name,
    within_saa_obsid,
    within_saa_ra,
    within_saa_dec,
    within_saa_start_time,
    within_saa_proc_time,
    within_saa_xrt_exp,
    within_saa_uvot_exp,
    within_saa_bat_exp,
    within_saa_archv_date]  # Each entry is a list itself
    for x in zip(*datalist1):
        outfile1.write('|{0}|{1}|{2}|{3}|{4}|{5}|{6}|{7}|{8}|{9}|\n'.format(*x))

with open(path+'swift_not_visible.txt', 'w') as outfile2:
    datalist2 = [not_visible_name,
    not_visible_obsid,
    not_visible_ra,
    not_visible_dec,
    not_visible_start_time,
    not_visible_proc_time,
    not_visible_xrt_exp,
    not_visible_uvot_exp,
    not_visible_bat_exp,
    not_visible_archv_date]  # Each entry is a list itself
    for x in zip(*datalist2):
        outfile2.write('|{0}|{1}|{2}|{3}|{4}|{5}|{6}|{7}|{8}|{9}|\n'.format(*x))

with open(path+'swift_visible.txt', 'w') as outfile3:
    datalist3 = [visible_name,
    visible_obsid,
    visible_ra,
    visible_dec,
    visible_start_time,
    visible_proc_time,
    visible_xrt_exp,
    visible_uvot_exp,
    visible_bat_exp,
    visible_archv_date]  # Each entry is a list itself
    for x in zip(*datalist3):
        outfile3.write('|{0}|{1}|{2}|{3}|{4}|{5}|{6}|{7}|{8}|{9}|\n'.format(*x))
