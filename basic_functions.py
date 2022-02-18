# define the right path to the data folders in this module
from config import *


# redshifts
# Create a dictionary with the redshift of each galaxy read from a special file
# File was created with the script "get_redshifts.py" in the folder /INPUT_DATA
# The name of the file is redshifts.txt and is located in /INPUT_DATA
# Its location was passed in config.py
redshifts = {}
with open(redshifts_file) as fh:
    for line in fh:
        gal, z = line.split()
        redshifts[gal] = z

# Read sky value
def read_sky(gal, band):
    if gal == "m3":
        gal = "3"
    elif gal == "m943":
        gal = "943"
    sky = np.loadtxt(folder_sky+"sky_"+gal+"_"+band+".txt")
    return sky


# Band '850' uses the suffix 'lp' instead of 'w'
def read_hst(gal,band,case='60'):
    if band in ['105','125','140','160','435'] and case=='30':
        print("WARNING: some bands only have data in 60 mas")
        print("Retrieving 60 mas data for gal,band = ", gal, ", ", band)
        case = '60'
    try:
         hdul = fits.open(folder_hst+gal+'_'+case+'mas_f'+band+"w_sci.fits")
    except:
        hdul = fits.open(folder_hst+gal+'_'+case+'mas_f'+band+"lp_sci.fits")
    return hdul


def read_muse_clean(gal):
    """Read the OII data cubes. Some data has a different path ('mos'). It also
       corrects one field in the header that is wrong in some cases.
    """

    z = redshifts[gal]
    combo = "udf10_c042_e031_"+gal+"_o2_Z_"+z
    if not os.path.isdir(folder_muse+'/'+combo):
        combo = 'udf_mos_c042_e030_'+gal+'_o2_Z_'+z
    hdul1 = fits.open(folder_muse+combo+"/"+combo+"_ssmooth_flux_common_OII3726_mclean5.0.fits")
    hdul2 = fits.open(folder_muse+combo+"/"+combo+"_ssmooth_flux_common_OII3729_mclean5.0.fits")
    hdul1[0].header["WCSAXES"] = 2
    hdul2[0].header["WCSAXES"] = 2
    return hdul1, hdul2


def total_lineflux(hdul1,hdul2,filename=None):
    """Sum of the flux in both lines of the doublet"""

    flux = hdul1[0].data + hdul2[0].data
    if filename:
        fits.writeto(filename,flux,hdul1[0].header,overwrite=True)
    return flux


# Some data is stored with a different filename structure
def read_snr_map(gal):
    z = redshifts[gal]
    try:
        combo = "udf10_c042_e031_"+gal+"_o2_Z_"+z
        hdul = fits.open(folder_muse+combo+"/"+combo+"_ssmooth_snr_common.fits")
    except:
        combo = "udf_mos_c042_e030_"+gal+"_o2_Z_"+z
        hdul = fits.open(folder_muse+combo+"/"+combo+"_ssmooth_snr_common.fits")
    return hdul


# WARNING: assuming FITS standard, pixels start at 1 (not at 0 !)
def get_wcs_coordinates(pixcrd, header):
    """ You provide a 2-column vector with the position of some pixels and you
        get the corresponding WCS coordinates of those points
    """

    w = wcs.WCS(header)
    world = w.wcs_pix2world(pixcrd, 1)
    return world


def get_pixel_coordinates(world,header):
    """ You provide a 2-column vector with the WCS position of some pixels and
        you get the corresponding pixel-coordinates of those points
    """

    w = wcs.WCS(header)
    pixcrd = w.wcs_world2pix(world, 1)
    return pixcrd


def oversampled_positions(image, oversampling):
    """
        Oversampled grid of pixels on top of an image by a given factor
        The positions of the small pixels are given in the units or frame of
        reference of the big pixels

        Warning: This is assuming that the center of the first big pixels is 1,
        maybe related to WCS transformations convention. 
        The left/bottom edges of the array are then in the position 0,5.
        The length of each small pixel is (1./oversampling), and we need to add
        just half of this quantity to reach the center of the first small pixel
        Then we continue adding integer amounts of the small pixel size

        Parameters
        ----------
        image : ndarray
            defines the original grid 
        oversampling : int
            factor of oversampling, the ratio of (large/small) pixels lengths

        Returns
        -------
        x : ndarray
            position of the small pixels in the reference of the large pixels
        y : ndarray
            position of the small pixels in the reference of the large pixels
        """
    
    npix = image.shape[0]
    x = 0.5 + (1./oversampling)/2 +  np.arange(0, npix * oversampling, 1)/float(oversampling)
    npix = image.shape[1]
    y = 0.5 + (1./oversampling)/2 +  np.arange(0, npix * oversampling, 1)/float(oversampling)
    
    return x, y


def native_pixel(band, case):
    """Size of the pixel in arcsec"""

    if band in ['105','125','140','160','435']:
        return 0.06
    elif band in ['606','775','814','850']:
        if case == '30':
            return 0.03
        elif case == '60':
            return 0.06


def direct_convolution(data,FWHM):
    """Convolution of 2D data with a gaussian given by its FWHM in pixels
    Arrays are enlarged by one FWHM in every direction before convolution,
    padding the data with zeros as necessary.
    """

    sigma = FWHM / (2 * np.sqrt(2 * np.log(2)))
    # FWHM must be passed in units of the pixels of the data
    FWHM = np.int(FWHM)+1
    data2 = np.zeros((data.shape[0] + 2 * FWHM, data.shape[1] + 2 * FWHM))
    data2[FWHM:data.shape[0] + FWHM, FWHM:data.shape[1] + FWHM] = data
    y, x = np.indices((data2.shape))
    # create gaussian PSF
    psf = np.exp(-((x - data2.shape[1] / 2) ** 2 + (y - data2.shape[0] / 2) ** 2) / (2.0 * sigma ** 2))
    psf /= psf.sum()  # normalisation PSF
    psf_shift = fftshift(psf)
    data_conv = irfft2(rfft2(data2) * rfft2(psf_shift))
    # get rid of the extra size added to the arrays
    data_conv = data_conv[FWHM:data.shape[0] + FWHM, FWHM:data.shape[1] + FWHM]
    return data_conv


def replace_negatives(im):
    mask = (im < 0)
    im[mask] = 0
    return im


def determine_smoothing(band1, band2, mode):
    """Calculate the 1-sigma size of a PSF to lower the resolution of images.
    * In mode 'equate' it will give 0 for the lowest resolution and the quadratic
    difference of both sigmas for the other band.
    *In mode 'smooth' it takes the high-resolution one to the resolution of the
    poor, and that one to the 0.6 arcsec resolution.
    * A return of 0 means no convolution is needed; be aware that trying to 
    convolve with a null gaussian is CATASTROPHIC
    """

    sigma1 = sigma_hst[band1]
    sigma2 = sigma_hst[band2]

    if mode == 'equate':
        # if both sigmas are the same there is no need to equate resolutions
        # Warning: convolving with null gaussian is catastrophic
        if (sigma1 == sigma2):
            psf1,psf2 = 0,0
        elif sigma1 > sigma2:
            psf1 = 0
            psf2 = np.sqrt(sigma1**2 - sigma2**2)
        else:
            psf2 = 0
            psf1 = np.sqrt(sigma2**2 - sigma1**2)

    elif mode == 'smooth':
        if (sigma1 == sigma2):
            print('Warning: method=smooth was supposed to be used with images \
                   of different resolutions, and these seem to have the same. \
                   They both will be lowered to 0.6 arcsec.')
            psf1 = np.sqrt(0.6**2 - sigma1**2)
            psf2 = np.sqrt(0.6**2 - sigma2**2)
        if sigma1 > sigma2:
            psf2 = np.sqrt(sigma1**2 - sigma2**2)
            psf1 = np.sqrt(0.6**2 - sigma1**2)
        else:
            psf1 = np.sqrt(sigma2**2 - sigma1**2)
            psf2 = np.sqrt(0.6**2 - sigma2**2)

    return psf1, psf2


def define_psf_muse_hst(gal, band1, band2):
    """To combine HST images and MUSE flux maps, here we calculate the extra
    smoothing needed as the quadratic difference between MUSE psf and the worst
    of both HST psfs
    """

    # These resolutions are defined in "config_paths.py", they are in arcsec
    sigma1 = sigma_hst[band1]
    sigma2 = sigma_hst[band2]

    # The pixel size for MUSE is 0.2 arcsec
    psf_m = psf_muse[gal] * 0.2 

    # An extra smoothing of 2 pixels is needed (ask Benoit)
    smooth = 2 * 0.2

    FWHM = np.sqrt(psf_m**2 + smooth**2)

    if (sigma1 >= sigma2):
        FWHM = np.sqrt(FWHM**2 - sigma1**2)
    else:
        FWHM = np.sqrt(FWHM**2 - sigma2**2)

    return FWHM


# undersamples an image to a coarser grid of pixels, considering a given factor
# The flux of a pixel in the new frame is the mean of the flux of all the small
# pixels that lie inside the area of the final big pixel
def rebin_image(data, oversampling):
    data_reshape = data.reshape(int(data.shape[0] / oversampling), oversampling, int(data.shape[1] / oversampling), oversampling)
    data_rebin = data_reshape.mean(1).mean(2)
    return data_rebin


def create_interpolant(x,y,values,fill_value=0):
    spl = scipy.interpolate.interp2d(x,y,values, fill_value=fill_value, kind='linear')
    return spl


###----- THIS IS THE MOST IMPORTANT FUNCTION OF THIS MODULE -----###
###-----                                                    -----###
def oversample_version(
                       hdul3, # expects a HST image
                       hdul1, # muse OII flux map
                       oversampling,
                       dx=0,
                       dy=0,
                       filename=None,
                       fill_value=0
    ):
    """It creates an oversampled version of the pixel grid of hdul1, and maps
       the image hdul3 to this refined grid by linear interpolation. This
       mapping is done in WCS coordinates, so the versions greated this way are
       physically aigned and they can be combined meaningfully.

       The shifts dx, dy are used in case we want to correct a bad astrometric
       calibration. They must be given in units of the hdul3 pixels.

       For the purpose of this job we normally expect hdul3 to be an HST image
       which is to be mapped to an oversampled grid of its corresponding OII
       flux map (hdul1).

       * "filename" is passed if you want to store the resulting image
       * "fill_value" is used to fill pixels that lie outside the original area,
       so that they can not be properly defined by the interpolation.
    """

    # position of the oversampled pixels in MUSE pixel units
    x,y = oversampled_positions(hdul1[0].data, oversampling)
    xx,yy = np.meshgrid(x,y)
    # 2-column vector, where each row represents a pixel of the oversampled grid
    pixcrd = np.column_stack((np.ravel(yy),np.ravel(xx)))
    # Lets get each pixel in this grid in WCS coordinates
    world_muse = get_wcs_coordinates(pixcrd, hdul1[0].header)
    # And now lets get them in coordinates of pixels in the hst frame of ref
    pixcrd2 = get_pixel_coordinates(world_muse, hdul3[0].header)
    # the interpolant is created directly in hst pixel coordinates to avoid artifacts
    npix = hdul3[0].data.shape[0]
    x_hst = np.arange(1,npix+1,1).astype(np.float)+dx
    npix = hdul3[0].data.shape[1]
    y_hst = np.arange(1,npix+1,1).astype(np.float)+dy
    spl = create_interpolant(x_hst,y_hst, hdul3[0].data, fill_value=fill_value)
    znew = np.zeros([x.size,y.size])
    # for some reason this was safer than trying directly a 2D interpolation 
    # It is quite fast anyways, so I advice to not spend time on this
    for i in range(y.size):
        for j in range(x.size):
            znew[j,i] = spl(pixcrd2[x.size*i+j,0],pixcrd2[x.size*i+j,1])
    # store the result
    if filename:
        hdr = oversampled_header(hdul1[0].header,oversampling)
        fits.writeto(filename,znew,hdr,overwrite=True)
    return znew


def oversampled_header(header,oversampling):
    # copy and edit the MUSE header
    hdr = header.copy()
    hdr["NAXIS1"] = header["NAXIS1"] * oversampling
    hdr["NAXIS2"] = header["NAXIS2"] * oversampling
    hdr["CD1_1"] = header["CD1_1"] / float(oversampling)
    hdr["CD2_2"] = header["CD2_2"] / float(oversampling)
    new_crpix1 = 0.5 + (header["CRPIX1"] - 0.5) * oversampling
    new_crpix2 = 0.5 + (header["CRPIX2"] - 0.5) * oversampling
    hdr["CRPIX1"] = new_crpix1
    hdr["CRPIX2"] = new_crpix2
    return hdr


def rebin_header(header,oversampling):
    # copy and edit the MUSE header
    hdr = header.copy()
    hdr["NAXIS1"] = header["NAXIS1"] / float(oversampling)
    hdr["NAXIS2"] = header["NAXIS2"] / float(oversampling)
    hdr["CD1_1"] = header["CD1_1"] * float(oversampling)
    hdr["CD2_2"] = header["CD2_2"] * float(oversampling)
    new_crpix1 = 0.5 + (header["CRPIX1"] - 0.5) / float(oversampling)
    new_crpix2 = 0.5 + (header["CRPIX2"] - 0.5) / float(oversampling)
    hdr["CRPIX1"] = new_crpix1
    hdr["CRPIX2"] = new_crpix2
    return hdr

