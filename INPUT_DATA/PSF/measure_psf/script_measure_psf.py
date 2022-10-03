import sys
sys.path.append("../")
from functions_psf import *

for star in ["1","2","3","4","5"]:
    for band in ["105", "125", "140", "160","435","606", "775", "814", "850"]:
        case = "60"
        do_all(star,band,case)
    for band in ["606", "775", "814", "850"]:
        case = "30"
        do_all(star,band,case)



