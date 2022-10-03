import os
import sys
import cv2
import pickle

import numpy as np
from astropy import wcs
import scipy.interpolate
from astropy.io import fits
from shutil import copyfile
from scipy import interpolate
import matplotlib.pyplot as plt
from scipy.optimize import least_squares
from numpy.fft import fftshift, rfft2, irfft2
from mpl_toolkits.axes_grid1 import make_axes_locatable


# define the right path to the data folders
folder_muse = "./INPUT_DATA/muse/"
folder_hst = "./INPUT_DATA/new_stamps/"
redshifts_file = "./INPUT_DATA/redshifts.txt"
folder_sky = "./INPUT_DATA/NEW_SKY/"
segmendir = "./INPUT_DATA/segment_maps/"


# Measured in pixels of MUSE data - (pixel = 0.2 arcsec)
psf_muse = {"1":3.09,"2":3.24,"3":3.16,"4":3.08,"5":3.25,"7":3.25,"8":2.99,
            "11":2.90,"12":2.92,"13":2.98,"14":3.08,"15":3.22,"16":2.99,
            "17":3.03,"21":2.99,"25":3.05,"26":2.97,"28":3.02,"30":2.99,"32":2.79,
            "37":2.94,"912":3.25,"919":2.95,"937":3.03,"943":3.19,"982":2.93,"1002":2.97}


# PSFs in arcsec measured in the 60mas images - FWHM
# I will assume it is the same for the 30 mas images... They likely are not !!!
# So this is something to be fixed/improved
sigma_hst =  {'105':0.208,
              '125':0.211,
              '140':0.219,
              '160':0.220,
              '435':0.110,
              '606':0.134,
              '775':0.123,
              '814':0.129,
              '850':0.120} # this numbers are reported in arcsec


