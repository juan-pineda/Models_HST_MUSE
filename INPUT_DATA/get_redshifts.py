import sys
import os
from astropy.io import fits

musedir = "./muse/"

# Get the name of all folders for muse data
directories = [x for x in os.listdir(musedir) if x[0] ==  "u"]

# Extract all pairs (galaxy-redshift) and write them to a file
file = open('redshifts.txt','w+')
for line in directories:
    # some directories are named, e.g, udf10_c042_e031_7382_o2_Z_1.093706
    if line[3]=="1":
        file.write(line.split('_')[3]+'\t'+line.split('_')[6]+'\n')
    # other directories use names like udf_mos_c042_e030_912_o2_Z_0.618915
    elif line[3]=="_":
        file.write(line.split('_')[4]+'\t'+line.split('_')[7]+'\n')
        
file.close()

