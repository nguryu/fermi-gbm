# This program takes the (RA, DEC) of an SGR recorded from the Fermi GBM Trigger catalogue and
# checks if it was visible to Fermi on that day.

import os
from gbm.time import Met
from gbm.data.poshist import PosHist

# Inputs
path = '/Users/RonnyNguyen/Desktop/FermiGBMSummer2020/data_dir/'
catalogue_file = path + 'fermi_sgr_list.txt'

# ================================================== #
# *** SGR LIST *** #
# ================================================== #
trigger_name = []
name = []
ra = []
dec = []
time = []
trigger_time = []
trigger_type = []
reliability = []

# Read in data.
with open(catalogue_file, 'r') as data_file:
    for line in data_file:
        line = line.strip()  # Remove whitespace.
        if not line:  # Skip empty lines.
            continue
        if not line.startswith('#'):
            col = line.split('|')
            trigger_name.append(col[1])
            name.append(col[2])
            ra.append(col[3])
            dec.append(col[4])
            time.append(col[5])
            trigger_time.append(col[6])
            trigger_type.append(col[7])
            reliability.append(col[8])

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
# Convert UTC time to mission elapsed time (MET in seconds) for Fermi.
met = []
for i in range(len(name)):
    time_conv = trigger_time[i][0:10] + 'T' + trigger_time[i][11:23]
    met.append(Met.from_iso(time_conv).met)

yy = trigger_time[0][2:4]
mm = trigger_time[0][5:7]
dd = trigger_time[0][8:10]

# Check if source is visible to Fermi or not.
within_saa_trigger_name = []
within_saa_name = []
within_saa_ra = []
within_saa_dec = []
within_saa_time = []
within_saa_trigger_time = []
within_saa_trigger_type = []
within_saa_reliability = []
within_saa_met = []
not_visible_trigger_name = []
not_visible_name = []
not_visible_ra = []
not_visible_dec = []
not_visible_time = []
not_visible_trigger_time = []
not_visible_trigger_type = []
not_visible_reliability = []
not_visible_met = []
visible_trigger_name = []
visible_name = []
visible_ra = []
visible_dec = []
visible_time = []
visible_trigger_time = []
visible_trigger_type = []
visible_reliability = []
visible_met = []

for i in range(len(name)):
    if trigger_time[i][0:10] >= '2009-06-01':  # Position files before mid-2009 would give KeyErrors.
        yy = trigger_time[i][2:4]
        mm = trigger_time[i][5:7]
        dd = trigger_time[i][8:10]

        # Check if v00, v01, v02 file.
        if os.path.isfile(path+'poshist_data/glg_poshist_all_'+yy+mm+dd+'_v00.fit') == True:
            poshist_file = path+'poshist_data/glg_poshist_all_'+yy+mm+dd+'_v00.fit'
        elif os.path.isfile(path+'poshist_data/glg_poshist_all_'+yy+mm+dd+'_v01.fit') == True:
            poshist_file = path+'poshist_data/glg_poshist_all_'+yy+mm+dd+'_v01.fit'
        else:
            poshist_file = path+'poshist_data/glg_poshist_all_'+yy+mm+dd+'_v02.fit'

        poshist = PosHist.open(poshist_file)

        # Check if inside SAA +/- 10 sec pad.
        # SAA = South Atlantic Anomaly where Fermi-GBM turns off due to the high particle flux.
        in_saa = poshist.get_saa_passage(met[i]) | poshist.get_saa_passage(met[i] + 10) | poshist.get_saa_passage(met[i] - 10)

        # Check if true location is occulted.
        occulted = poshist.location_visible(ra_deg[i], dec_deg[i], met[i]) == False

        if in_saa:
          # print("Fermi GBM is in SAA +/- 10 sec")
          within_saa_trigger_name.append(trigger_name[i])
          within_saa_name.append(name[i])
          within_saa_ra.append(ra[i])
          within_saa_dec.append(dec[i])
          within_saa_time.append(time[i])
          within_saa_trigger_time.append(trigger_time[i])
          within_saa_trigger_type.append(trigger_type[i])
          within_saa_reliability.append(reliability[i])
          within_saa_met.append(met[i])
        elif occulted:
          # print("Position not visible to Fermi GBM")
          not_visible_trigger_name.append(trigger_name[i])
          not_visible_name.append(name[i])
          not_visible_ra.append(ra[i])
          not_visible_dec.append(dec[i])
          not_visible_time.append(time[i])
          not_visible_trigger_time.append(trigger_time[i])
          not_visible_trigger_type.append(trigger_type[i])
          not_visible_reliability.append(reliability[i])
          not_visible_met.append(met[i])
        else:
          # print("Source should be visible to Fermi GBM")
          visible_trigger_name.append(trigger_name[i])
          visible_name.append(name[i])
          visible_ra.append(ra[i])
          visible_dec.append(dec[i])
          visible_time.append(time[i])
          visible_trigger_time.append(trigger_time[i])
          visible_trigger_type.append(trigger_type[i])
          visible_reliability.append(reliability[i])
          visible_met.append(met[i])

# Write relevant data into new file
with open(path+'fermi_within_saa.txt', 'w') as outfile1:
    datalist1 = [within_saa_trigger_name,
    within_saa_name,
    within_saa_ra,
    within_saa_dec,
    within_saa_time,
    within_saa_trigger_time,
    within_saa_trigger_type,
    within_saa_reliability,
    within_saa_met]  # Each entry is a list itself
    outfile1.write('# col 0-1 = trigger_name, name\n')
    outfile1.write('# col 2-3 = ra, dec\n')
    outfile1.write('# col 4-5 = time, trigger_time\n')
    outfile1.write('# col 6 = trigger_type\n')
    outfile1.write('# col 7 = reliability\n')
    outfile1.write('# col 8 = met\n')
    for x in zip(*datalist1):
        outfile1.write('|{0}|{1}|{2}|{3}|{4}|{5}|{6}|{7}|{8}\n'.format(*x))

with open(path+'fermi_not_visible.txt', 'w') as outfile2:
    datalist2 = [not_visible_trigger_name,
    not_visible_name,
    not_visible_ra,
    not_visible_dec,
    not_visible_time,
    not_visible_trigger_time,
    not_visible_trigger_type,
    not_visible_reliability,
    not_visible_met]
    outfile2.write('# col 0-1 = trigger_name, name\n')
    outfile2.write('# col 2-3 = ra, dec\n')
    outfile2.write('# col 4-5 = time, trigger_time\n')
    outfile2.write('# col 6 = trigger_type\n')
    outfile2.write('# col 7 = reliability\n')
    outfile2.write('# col 8 = met\n')
    for x in zip(*datalist2):
        outfile2.write('|{0}|{1}|{2}|{3}|{4}|{5}|{6}|{7}|{8}\n'.format(*x))

with open(path+'fermi_visible.txt', 'w') as outfile3:
    datalist3 = [visible_trigger_name,
    visible_name,
    visible_ra,
    visible_dec,
    visible_time,
    visible_trigger_time,
    visible_trigger_type,
    visible_reliability,
    visible_met]
    outfile3.write('# col 0-1 = trigger_name, name\n')
    outfile3.write('# col 2-3 = ra, dec\n')
    outfile3.write('# col 4-5 = time, trigger_time\n')
    outfile3.write('# col 6 = trigger_type\n')
    outfile3.write('# col 7 = reliability\n')
    outfile3.write('# col 8 = met\n')
    for x in zip(*datalist3):
        outfile3.write('|{0}|{1}|{2}|{3}|{4}|{5}|{6}|{7}|{8}\n'.format(*x))
