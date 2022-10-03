# Models\_HST\_MUSE
Library for mocking MUSE OII line flux maps using HST images

Here are the most important functions to allow the creation of such models.
Please run the built in example to understand what are the main steps.
Here we provide sample data for running tests on 9 galaxies. Yet remember
that github is not intended for storing large datasets, so, for large samples
simply update the paths inside config.py to the folders with the muse
datacubes, HST images, etc.  

Here it is supposed that you have precomputed values of sky and psf for hubble
bands (see the sections below). We provide the necessary inputs for the sample
data. This repository still needs to be updated with some techniques to perform
a meaningful model fitting using statistical indicators. Meanwhile, the example
creates just a toy model to allow you grasping the most important features of
the code. 

# Measurements of the SKY

The procedure for measuring sky is now included in this repo. Check the script
"sky\_authomatic.py". Some functions specifically dedicated for this task were
added to the module "basic\_functions.py", clearly separated by a distinct mark
and title, below the body of functions already present in there.

# PSF determination

Take a look at the folder "INPUT\_DATA/PSF/", everything is self-contained there.




