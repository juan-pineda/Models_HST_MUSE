# With this script I open the PSF images one by one, in order to create a mask
# for each one. That is necessary so the fit to the gaussian later on is 
# performed only on significant pixels

star=$1

for band in 606 775 814
do
ds9 -log -zoom 3 -regions shape circle ../psf/${star}_30mas_f${band}w_sci.fits 
done
band=850
ds9 -log -zoom 3 -regions shape circle ../psf/${star}_30mas_f${band}lp_sci.fits

for band in 606 775 814 105 125 140 160 435
do
ds9 -log -zoom 3 -regions shape circle ../psf/${star}_60mas_f${band}w_sci.fits
done
band=850
ds9 -log -zoom 3 -regions shape circle ../psf/${star}_60mas_f${band}lp_sci.fits



