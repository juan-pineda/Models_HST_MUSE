

# These functions are brought from `TEST_3_masking_not_psf/test_PSF`
# from files: useful_psf_functions.py & preliminar_tests/test_1.py

from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
import pyregion
import astropy.convolution
#from funciones import *
from scipy.optimize import leastsq
from mpl_toolkits.axes_grid1 import make_axes_locatable


maskdir = "/home/juan/Desktop/Models_HST_MUSE/INPUT_DATA/PSF/masks"
psfdir = "/home/juan/Desktop/Models_HST_MUSE/INPUT_DATA/PSF/psf"
measuredir = "/home/juan/Desktop/Models_HST_MUSE/INPUT_DATA/PSF/measure_psf"

# Esta función carga las regiones definidas via ds9
# Y almacena las máscaras 2D, con 1s y 0s
def create_masks(star,band,case):
    r = pyregion.open(maskdir+"/mask_psf_"+star+"_"+band+"_"+case+"mas.reg")
    f = read_data_n_header(star,band,case)
    mask = r.get_mask(hdu=f)
    filename = maskdir+"/mask_psf_"+star+"_"+band+"_"+case+"mas.fits"
    fits.writeto(filename,mask*1,f.header,overwrite=True)

# case : pixel scale in arcsec, i.e., 30 or 60
def read_mask(star,band,case):
    filename = maskdir+"/mask_psf_"+star+"_"+band+"_"+case+"mas.fits"
    f = fits.open(filename)
    return f[0].data

# read the psf data
def read_data_n_header(star,band,case):
    if band == "850":
        f = fits.open(psfdir+"/"+star+"_"+case+"mas_f"+band+"lp_sci.fits")
    else:
        f = fits.open(psfdir+"/"+star+"_"+case+"mas_f"+band+"w_sci.fits")
    return f[0]


# read the psf data
def read_data(star,band,case):
    f = read_data_n_header(star,band,case)
    return f.data

def gaussian2D(tpl,image):
    """Returns a gaussian function with the given parameters"""
    gauss = np.zeros(image.shape)
    x = np.arange(image.shape[0])
    y = np.arange(image.shape[1])
    yy,xx = np.meshgrid(y,x)

    x0 = tpl[0]
    y0 = tpl[1]
    A = tpl[2]
    std = tpl[3]
    cc = tpl[4]
    gauss = cc + A*np.exp(-((xx-x0)**2+(yy-y0)**2)/(std**2)/2)
    return gauss

def init_guess_gauss2D(image):
    x0 = image.shape[0]/2
    y0 = image.shape[1]/2
    A = np.max(image)
    std = 3
    cc1 = np.mean(image[0,:])
    cc2 = np.mean(image[:,0])
    cc3 = np.mean(image[-1,:])
    cc4 = np.mean(image[:,-1])
    cc = (cc1+cc2+cc3+cc4)/4
    p0 = [x0,y0,A,std,cc]
    return p0


def ErrorFunc2(tpl,image,mask):
    model = gaussian2D(tpl,image)*mask
    data = image*mask
    chi = np.ravel(model-data)
    return chi

def do_all(star,band,case):
    data = read_data(star,band,case)
    mask = read_mask(star,band,case)
    model,std = gaussian_fit(data,mask)
    show_result(data*mask,model*mask,log=False)
    plt.savefig(measuredir+"/"+star+"_"+band+"_"+case+"mas_linear.jpg")
    show_result(data*mask,model*mask,log=True)
    plt.savefig(measuredir+"/"+star+"_"+band+"_"+case+"mas_log.jpg")
    print(star,"\t",band,"\t",case,"\t",std)


def show_result(data,model,log=False):

    mask = (data > 0)
    data[~mask] = np.NaN
    model[~mask] = np.NaN

    plt.close()
    fig = plt.figure(figsize=(12,4))
    ax0 = fig.add_subplot(131)
    ax1 = fig.add_subplot(132)
    ax2 = fig.add_subplot(133)
    ax0.set_axis_off()
    ax1.set_axis_off()
    ax2.set_axis_off()

    if log:
        data = np.log10(data)
        model = np.log10(model)

    vmin = np.percentile(data[~np.isnan(data)],10)
    vmax = np.percentile(data[~np.isnan(data)],99)

    resi = (data-model)
    resi_p = np.abs(100 * resi/data)

    im = ax0.imshow(data,vmin=vmin,vmax=vmax)
    divider = make_axes_locatable(ax0)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    fig.colorbar(im, ax=ax0,cax=cax)
    im = ax1.imshow(model,vmin=vmin,vmax=vmax)
    divider = make_axes_locatable(ax1)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    fig.colorbar(im, ax=ax1,cax=cax)
    im = ax2.imshow(resi_p,vmin=0,vmax=100)
    divider = make_axes_locatable(ax2)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    fig.colorbar(im, ax=ax2,cax=cax)

def gaussian_fit(image,mask):
    p0 = init_guess_gauss2D(image)
    best_params,success = leastsq(ErrorFunc2, p0, args = (image,mask))
    g1 = gaussian2D(best_params,image)
    std = best_params[3]
    return g1,std


def reduced_X2(data,model,Npars,std):
    # this is to ensure that minimum error is std dev
    model[(model < 0)] = 0
    sigma = np.sqrt(data + std**2)
    X2 = (np.ravel((data-model)/sigma)**2).sum()
    X2 = X2 / (data.size - Npars)
    return X2




