# Given a list of SGRs, this program extracts position history data from the
# Fermi database: https://heasarc.gsfc.nasa.gov/FTP/fermi/.

import os

# Inputs
path = '/Users/RonnyNguyen/Desktop/FermiGBMSummer2020/data_dir/'
catalogue_file = path + 'sgr_list.txt'

# Initialize lists.
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

date = []
for i in range(len(start_time)):
    if start_time[i][0:10] >= '2008-08-16':  # No Fermi data before August 16, 2008.
        date.append(start_time[i][0:10])

date = list(sorted(dict.fromkeys(date)))  # Sort dates in ascending order and remove duplicates.

# Download the position history data for the day of observation.
for i in range(len(date)):
        yy = date[i][2:4]
        mm = date[i][5:7]
        dd = date[i][8:10]
        link = "https://heasarc.gsfc.nasa.gov/FTP/fermi/data/gbm/daily/20"+yy+"/"+mm+"/"+dd+"/current/glg_poshist_all_"+yy+mm+dd+"_v00.fit"

        dest = path + "glg_poshist_all_"+yy+mm+dd+"_v00.fit"

        os.system("curl  %s > %s" % (link, dest))
