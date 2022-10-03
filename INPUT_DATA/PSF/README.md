
# PSF determination

So, we have sample HST images for five stars in the different bands and pixel
scales of the HST. The idea is simply fitting gaussian functions, but first we
mask the stellar footprints in order to make the fitting only over significant
pixels. This is necessary because the actual form of the PSF footprint is not
really gaussian, so, when including the whole area the outskirts of the image 
may bias the gaussian fits and produce weird values of the PSF.

**NOTE:** You must specify the specific path to the three subfolders in "./PSF/"
inside file "functions\_psf.py"

**IMPORTANT:** You need to install pyregion runnig this code snippet:  
\$ conda install -c conda-forge pyregion

Step by step:

- Folder "./psf/" contains the data of the stars that are used to measure the PSF

- The masks are created in the folder "./masks/". For that, you must first move 
inside that directory and run the script "launch\_ds9.sh", which also needs a
number between 1 and 5 as a paremeter, indicating one of the stars. This script
will display in ds9 all the images available for that star, one by one, and your
task is to draw a region enclosing the pixels you want to use for the gaussian
fit. Then this regions must be exported in .reg format using "image" as your
coordinate system, and "ds9" as your "file format". The filenames should follow
this pattern: "mask\_psf\_1\_105\_60mas.reg"

- Next step is to run the file "script\_mk\_masks.py" while staying still inside
the same folder (masks). This scripts will read the .reg files and convert the
information into .fits images of the masks.

- Now that the masks are ready you should move to the folder "./measure\_psf/".
There you first run the file "script\_measure\_psf.py". That will perform the
gaussian fit over the accepted area for each one of the cases. The resulting
sigma of all the psfs are prnted on screen and some diagnostic plots witht the
results of the gaussian modeling are stored in the same folder (measure\_psf).
All the results printed on screen are manually prompted to a text file named
"psf\_measurements.txt", with 4 columns: [star, band, pixel scale, sigma of psf]

- Finally run the file "script\_psf\_compilation.py". This script will report
the sigma of the PSF in each band after averaging over all the stars, and it will
also report an error as the standard deviation between them. These numbers are
printed on screen and we manually prompted them tot he file "psf\_compilation.txt"
inside the same folder (measure\_psf).











