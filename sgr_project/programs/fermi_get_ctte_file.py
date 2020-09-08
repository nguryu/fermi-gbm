# Given a list of SGRs, this program extracts TTE data from the
# Fermi database: https://heasarc.gsfc.nasa.gov/FTP/fermi/.
# Edit the date in line 42 to a specific date if you want to download TTE data for that one day.
# Edit the '==' operator to '>=' in line 42 if you want to download data continuously after the designated date.

import os

# Inputs
path = '/Users/RonnyNguyen/Desktop/FermiGBMSummer2020/data_dir/'
input_file = path + 'fermi_sgr_visible.txt'
output_path = '/Users/RonnyNguyen/Desktop/FermiGBMSummer2020/data_dir/ctte_data/'

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
    if trigger_time[i][0:10] == '2013-01-13':  # Non-standardized TTE data before 2013-01-01
        date.append(trigger_time[i][0:10])

date = list(sorted(dict.fromkeys(date)))  # Sort dates in ascending order and remove duplicates.

suff = ['b0', 'b1', 'n0', 'n1', 'n2', 'n3', 'n4', 'n5', 'n6', 'n7', 'n8', 'n9', 'na', 'nb']

hour = ['00', '01', '02', '03', '04', '05', '06', '07', '08', '09', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', '23']

# Download the TTE data for the day of observation.
for i in range(len(date)):
        yy = date[i][2:4]
        mm = date[i][5:7]
        dd = date[i][8:10]

        for j in range(len(suff)):
            for k in range(len(hour)):
                link = "https://heasarc.gsfc.nasa.gov/FTP/fermi/data/gbm/daily/20"+yy+"/"+mm+"/"+dd+"/current/glg_tte_"+suff[j]+"_"+yy+mm+dd+"_"+hour[k]+"z_v00.fit.gz"

                # Make file in directory for date.
                dir = os.path.join(output_path, "20"+yy+"-"+mm+"-"+dd)
                if not os.path.exists(dir):
                    os.mkdir(dir)

                dest = output_path + "/" + "20"+yy+"-"+mm+"-"+dd + "/" + "glg_tte_"+suff[j]+"_"+yy+mm+dd+"_"+hour[k]+"z_v00.fit.gz"

                os.system("curl %s > %s" % (link, dest))
