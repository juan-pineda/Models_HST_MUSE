import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt('psf_measurements.txt')

case = 60
for band in [105, 125, 140, 160, 435, 606, 775, 814, 850]:
    index = (data[:,2] == case) & (data[:,1] == band)
    psf = np.mean(data[:,3][index])
    err = np.std(data[:,3][index])
    print('PSF for case '+str(case)+' in band '+str(band)+' is '+str(np.round(psf,2))+' ('+str(np.round(err,2))+')  in pixels')

print('\n')

case = 30
for band in [606, 775, 814, 850]:
    index = (data[:,2] == case) & (data[:,1] == band)
    psf = np.mean(data[:,3][index])
    err = np.std(data[:,3][index])
    print('PSF for case '+str(case)+' in band '+str(band)+' is '+str(np.round(psf,2))+' ('+str(np.round(err,2))+')  in pixels')



