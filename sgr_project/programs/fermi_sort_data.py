# Given a list of SGRs, this program sorts them in ascending order by date of observation using the Pandas library.
# It also converts UTC time to MET time to make it more convenient to run the targeted search.
# Edit line 52 to sort by desired category.

import pandas as pd

# Inputs.
path = '/Users/RonnyNguyen/Desktop/FermiGBMSummer2020/data_dir/'
input_file = path + 'fermi_visible.txt'

# Initialize lists.
trigger_name = []
name = []
ra = []
dec = []
time = []
trigger_time = []
trigger_type = []
reliability = []
met = []

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
            met.append(col[9])

# Use Pandas module to put data into a Data Frame.
df = pd.DataFrame({'trigger_name':trigger_name,
                    'name':name,
                    'ra':ra,
                    'dec':dec,
                    'time':time,
                    'trigger_time':trigger_time,
                    'trigger_type':trigger_type,
                    'reliability':reliability,
                    'met':met})

# Sort columns by desired column, this creates a new Data Frame.
sorted = df.sort_values(by='trigger_time', ascending=False)

trigger_name = sorted['trigger_name'].tolist()
name = sorted['name'].tolist()
ra = sorted['ra'].tolist()
dec = sorted['dec'].tolist()
time = sorted['time'].tolist()
trigger_time = sorted['trigger_time'].tolist()
trigger_type = sorted['trigger_type'].tolist()
reliability = sorted['reliability'].tolist()
met = sorted['met'].tolist()

# Write relevant data into new file
with open(path+'fermi_sgr_sorted.txt', 'w') as outfile:
    datalist = [trigger_name, name, ra, dec, time, trigger_time, trigger_type, reliability, met]  # Each entry is a list itself.
    outfile.write('# col 0-1 = trigger_name, name\n')
    outfile.write('# col 2-3 = ra, dec\n')
    outfile.write('# col 4-5 = time, trigger_time\n')
    outfile.write('# col 6 = trigger_type\n')
    outfile.write('# col 7 = reliability\n')
    outfile.write('# col 8 = met\n')
    for x in zip(*datalist):
        outfile.write('|{0}|{1}|{2}|{3}|{4}|{5}|{6}|{7}|{8}|\n'.format(*x))
