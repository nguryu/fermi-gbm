# This program reads in the Swift Master Catalog and searches for SGR sources.
# The SGR sources, along with relevant information, are wirrten into a file.

# Inputs
path = '/Users/RonnyNguyen/Desktop/FermiGBMSummer2020/data_dir/'
catalogue_file = path + 'SwiftMasterCatalog.txt'
outfile_name = path + 'sgr_list.txt'

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

# Search for SGR sources.
sgr_name = []
sgr_obsid = []
sgr_ra = []
sgr_dec = []
sgr_start_time = []
sgr_proc_time = []
sgr_xrt_exp = []
sgr_uvot_exp = []
sgr_bat_exp = []
sgr_archv_date = []
for i in range(len(name)):
    if name[i][0:3] == 'SGR' and name[i][0:4] != 'SGRA':  # Exclude observations from Sagittarius A*.
        sgr_name.append(name[i])
        sgr_obsid.append(obsid[i])
        sgr_ra.append(ra[i])
        sgr_dec.append(dec[i])
        sgr_start_time.append(start_time[i])
        sgr_proc_time.append(proc_time[i])
        sgr_xrt_exp.append(xrt_exp[i])
        sgr_uvot_exp.append(uvot_exp[i])
        sgr_bat_exp.append(bat_exp[i])
        sgr_archv_date.append(archv_date[i])

# Write relevant data into new file
with open(outfile_name, 'w') as outfile:
    figlist = [sgr_name,
    sgr_obsid,
    sgr_ra,
    sgr_dec,
    sgr_start_time,
    sgr_proc_time,
    sgr_xrt_exp,
    sgr_uvot_exp,
    sgr_bat_exp,
    sgr_archv_date]  # Each entry is a list itself
    for x in zip(*figlist):
        outfile.write('|{0}|{1}|{2}|{3}|{4}|{5}|{6}|{7}|{8}|{9}|\n'.format(*x))
