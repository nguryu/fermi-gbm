# Given a list of SGRs, this program extracts CTTE data from the
# Fermi database: https://heasarc.gsfc.nasa.gov/FTP/fermi/.

import os

# Inputs
path = '/Users/RonnyNguyen/Desktop/FermiGBMSummer2020/data_dir/'
input_file = path + 'visible.txt'
output_file = '/Users/RonnyNguyen/Desktop/FermiGBMSummer2020/data_dir/ctte_data/'

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
with open(input_file, 'r') as data_file:
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
    if start_time[i][0:10] == '2020-06-19':  # No CTTE data before August 10, 2010.
        date.append(start_time[i][0:10])

date = list(sorted(dict.fromkeys(date)))  # Sort dates in ascending order and remove duplicates.

suff = ['b0', 'b1', 'n0', 'n1', 'n2', 'n3', 'n4', 'n5', 'n6', 'n7', 'n8', 'n9', 'na', 'nb']

num = ['00', '01', '02', '03', '04', '05', '06', '07', '08', '09', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', '23']

# Download the position history data for the day of observation.
for i in range(len(date)):
        yy = date[i][2:4]
        mm = date[i][5:7]
        dd = date[i][8:10]

        for j in range(len(suff)):
            for k in range(len(num)):
                link = "https://heasarc.gsfc.nasa.gov/FTP/fermi/data/gbm/daily/20"+yy+"/"+mm+"/"+dd+"/current/glg_tte_"+suff[j]+"_"+yy+mm+dd+"_"+num[k]+"z_v00.fit.gz"

                # Make file in directory for date.
                dir = os.path.join(output_file, yy+mm+dd)
                if not os.path.exists(dir):
                    os.mkdir(dir)

                dest = output_file + "/" + yy+mm+dd + "/glg_tte_"+suff[j]+"_"+yy+mm+dd+"_"+num[k]+"z_v00.fit.gz"

                os.system("curl  %s > %s" % (link, dest))
