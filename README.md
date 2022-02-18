# Models_HST_MUSE
Library for mocking MUSE OII line flux maps using HST images

Here are the most important functions to allow the creation of such models.
Please run the built in example to understand what are the main steps.
Here we provide sample data for running tests on 9 galaxies. Yet remember
that github is not intended for storing large datasets, so, for large samples
simply update the paths inside config.py to the folders with the muse
datacubes, HST images, etc.  

Here it is supposed that you have precomputed values of sky and psf for hubble
bands. We provide the necessary inputs for the sample data. In a couple of days
this repository will be updated with extra features and explanations for doing
that, and also to perform the true model fitting with serious statistical
meaning. Meanwhile, the example creates just a toy model to allow you grasping
the most important features of the code, so later on you can go as far as you
want.  


