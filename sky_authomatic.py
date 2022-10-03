from basic_functions import *

galaxies = ["3","15","37","912","919","937","943","982","1002"]
bands = ["105","125","140","160","435","606","775","814","850"]

for gal in galaxies:
    print("working in gal ",gal)

    # read in the data cubes
    hdul1, hdul2 = read_muse_clean(gal)

    for band in bands:

        # read in the hubble images
        hdul3 = read_hst(gal,band,case='60')

        # get the positions of the  SKY pixels only
        mask, header, seg_mask = mask_sky(gal)
        new_mask = interpolate_mask(mask,header,hdul3)

        if band in ["435","606","775","814","850"]:
            # this filename is to store a histogram of the sky
            filename = folder_sky+"histogram_"+gal+"_"+band+".png"
            sky, std = get_sky_n_fluctuations(hdul3, new_mask, graph=True, filename=filename)

        # for these bands there are some border effects
        elif band in ["105","125","140","160"]:
            erosion = erode_sky(new_mask,n=5,iterations=5)
            filename = folder_sky+"histogram_"+gal+"_"+band+".png"
            sky, std = get_sky_n_fluctuations(hdul3, erosion, graph=True, filename=filename)

        filename = folder_sky+"sky_"+gal+"_"+band+".txt"
        np.savetxt(filename,np.array([sky]))

        hdul3.close()
    hdul1.close()
    hdul2.close()

