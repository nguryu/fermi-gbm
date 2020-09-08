# Given a list of SGRs, this program extracts position history data from the
# Fermi database: https://heasarc.gsfc.nasa.gov/FTP/fermi/.
# Edit the date in line 41 to a specific date if you want to download TTE data for that one day.
# Edit the '==' operator to '>=' in line 41 if you want to download data continuously after the designated date.

import os

# Inputs
path = '/Users/RonnyNguyen/Desktop/FermiGBMSummer2020/data_dir/'
input_file = path + 'fermi_sgr_list.txt'

# Initialize lists.
trigger_name = []
name = []
ra = []
dec = []
time = []
trigger_time = []
trigger_type = []
reliability = []

# Read in data.
with open(input_file, 'r') as data_file:
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

date = []
for i in range(len(trigger_time)):
    if trigger_time[i][0:10] >= '2008-08-16':  # No Fermi data before August 16, 2008.
        date.append(trigger_time[i][0:10])

date = list(sorted(dict.fromkeys(date)))  # Sort dates in ascending order and remove duplicates.

# Download the position history data for the day of observation.
for i in range(len(date)):
        yy = date[i][2:4]
        mm = date[i][5:7]
        dd = date[i][8:10]

        poshist_file_v00 = path+'poshist_data/glg_poshist_all_'+yy+mm+dd+'_v00.fit'
        poshist_file_v01 = path+'poshist_data/glg_poshist_all_'+yy+mm+dd+'_v01.fit'

        # Check if v00 or v01 file already exists in local directory.
        if os.path.isfile(poshist_file_v00) == True or os.path.isfile(poshist_file_v01) == True:
            continue
        elif os.path.isfile(poshist_file_v00) == False and os.path.isfile(poshist_file_v01) == False:
            link = "https://heasarc.gsfc.nasa.gov/FTP/fermi/data/gbm/daily/20"+yy+"/"+mm+"/"+dd+"/current/glg_poshist_all_"+yy+mm+dd+"_v00.fit"
            dest = path + "poshist_data/glg_poshist_all_"+yy+mm+dd+"_v00.fit"

        os.system("curl  %s > %s" % (link, dest))
