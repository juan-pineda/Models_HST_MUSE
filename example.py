from config import *
import basic_functions as bf

# Data for running examples with these galaxies is available in /INPUT_DATA/
galaxies = ["3","15","37","912","919","937","943","982","1002"]

# In this example we will create a toy model combinning 2 bands of HST to mock
# the corresponding OII flux map from MUSE.

# These are the steps:

# - OII total flux map is created summing the flux in both lines
# - HST images are sky subtracted and equated in resolution
# - Each image is interpolated to a grid that is muse-like but oversampled
# - A linear combination is done with fixed weights; one positive, one negative
# - Negative values are replaced by zeroes
# - The image is smoothed to muse resolution, and rebinned to the same shape
# - Image is rescaled so the total flux is the same as the total OII flux 
# - Residual is calculated

# The byproducts are stored to allow further inspection, so they can help you
# to check the effect of every step. The methodology of model fitting is 
# intentionally naive, just to illustrate the basic capabilities of the code,
# on top of which you can elaborate further once you get familiar with it.


###--------------------------------------###
###--- BEGIN OF MANUAL CONFIGURATIONS ---###

# Let's create some example products for galaxy 3
# Choose one band 'blue' and one band 'red'
gal = "3"
band1 = "160"
band2 = "814"
# Weights for the linear combination of both bands
w1 = 1
w2 = -0.2

# We will create a grid on top of the field-of-view covered by MUSE
# but oversampling each pixel by a factor of 8
oversampling = 8
# Folder to store the results
outdir = './results_example'
# Let's work with the HST images with pixels of 30 mas
case = '30'

###--- END OF MANUAL CONFIGURATIONS ---###
###------------------------------------###

print("Example with the next configuration:")
print("galaxy: ", gal)
print("Bands: ", band1, ", ", band2)
print("Weights: ", w1, ", ", w2)
print("Working with images of ", case, "miliarcsec pixel size (mas)")
print("\n ... \n")

###------------------------###
###--- Read the data in ---###

print("Reading data...")

# Read the 2 flux maps in the OII lines of the doublet
hdul1, hdul2 = bf.read_muse_clean(gal)

# Sums both maps to have the total flux of the doublet. Store the result.
filename = os.path.join(outdir,"gal_"+gal+"_data.fits")
data = bf.total_lineflux(hdul1, hdul2, filename)
print("Total flux map stored in: ", filename)

# Read the 2 HST images and subtract sky
hdul3 = bf.read_hst(gal, band1, case)
h3 = hdul3[0].data - bf.read_sky(gal, band1)
hdul4 = bf.read_hst(gal, band2, case)
h4 = hdul4[0].data - bf.read_sky(gal, band2)

###-----------------------------------------------###
###--- Equate the resolutions and subtract sky ---###

# check the necessary smoothing to equate resolutions
psf1, psf2 = bf.determine_smoothing(band1, band2, mode='equate')

if psf1 > 0:
    # fwhm resolution is divided by pixel-size to have it in numer of pixels
    psf1 = psf1 / bf.native_pixel(band1, case)
    h3 = bf.direct_convolution(h3, psf1)

if psf2 > 0:
    psf2 = psf2 / bf.native_pixel(band2, case)
    h4 = bf.direct_convolution(h4, psf2)

filename = os.path.join(outdir,"gal_"+gal+"_hst_"+band1+"_input.fits")
fits.writeto(filename, h3, hdul3[0].header, overwrite=True)
print("Input image stored in: ", filename)
filename = os.path.join(outdir,"gal_"+gal+"_hst_"+band2+"_input.fits")
fits.writeto(filename, h4, hdul4[0].header, overwrite=True)
print("Input image stored in: ", filename)

###---------------------------------------------------------------###
###--- Interpolate the HST images to the oversampled MUSE grid ---###

# dx, dy are for small shifts in case of astrometry errors. We set them to 0
# fill_value is for pixels outside original area, where interpolation is not valid
filename = os.path.join(outdir,"gal_"+gal+"_hst_"+band1+"_interpolated.fits")
h3 = bf.oversample_version(
                        hdul3, 
                        hdul1, 
                        oversampling,
                        dx=0,
                        dy=0,
                        filename=filename,
                        fill_value=0
    )
print("Oversampled image stored in: ", filename)

filename = os.path.join(outdir,"gal_"+gal+"_hst_"+band2+"_interpolated.fits")
h4 = bf.oversample_version(
                        hdul4, 
                        hdul1, 
                        oversampling,
                        dx=0,
                        dy=0,
                        filename=filename,
                        fill_value=0
    )
print("Oversampled image stored in: ", filename)

###---------------------------------------------------------------###
###--- Linear combination and reduction for comparing to muse ----###

model = w1*h3 + w2*h4
model = bf.replace_negatives(model)

filename = os.path.join(outdir,"gal_"+gal+"_highres_model.fits")
hdr = bf.oversampled_header(hdul1[0].header,oversampling)
fits.writeto(filename, model, hdr, overwrite=True)
print("High resolution model stored in: ", filename)

# smooth to muse resolution
psf = bf.define_psf_muse_hst(gal, band1, band2)
model = bf.direct_convolution(model, psf)

filename = os.path.join(outdir,"gal_"+gal+"_smoothed_model.fits")
fits.writeto(filename, model, hdr, overwrite=True)
print("Smoothed high-sampling model stored in: ", filename)

# resample
model = bf.rebin_image(model, oversampling)

filename = os.path.join(outdir,"gal_"+gal+"_rebinned_model.fits")
fits.writeto(filename, model, hdul1[0].header, overwrite=True)
print("Rebinned model stored in: ", filename)

#renormalization
model = model * data.sum() / model.sum()

filename = os.path.join(outdir,"gal_"+gal+"_model.fits")
fits.writeto(filename, model, hdul1[0].header, overwrite=True)
print("Input image stored in: ", filename)

residual = data - model

filename = os.path.join(outdir,"gal_"+gal+"_residual.fits")
fits.writeto(filename, residual, hdul1[0].header, overwrite=True)
print("Input image stored in: ", filename)










 

