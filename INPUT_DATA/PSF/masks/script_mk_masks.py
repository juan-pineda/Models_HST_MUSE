import sys
sys.path.append("../")
from functions_psf import *

for star in ["1","2","3","4","5"]:
    for band in ["606", "775", "814", "850", "105", "125", "140", "160", "435"]:
        create_masks(star,band,case="60")
    for band in ["606", "775", "814", "850"]:
        create_masks(star,band,case="30")

