# Optical_Detectors_Assignment
This repository is my optical detectors assignment.

When the script optical_detectors.py is run, you can add arguments for the path to folders F336W and F555W (in that order). It should run fine if you do not, because
if paths are not specified the script should assume they are in the "data" folder of this github, and they are. I just personally think it makes more sense to specify your file path each time. 

# What the script will do

Using data from the F336W and F555W filter images, this script will median-combine each set of filter images to one per filter. This will get rid of (most) any cosmic rays. 
Then peaks are detected in one image and catagorized as stars or not (passing or failing checks). The list of peak coordinates is checked in the second image, then a final catalog of stars is produced of sources that look like stars at the same coordinates in both images.

Photometery is performed to get the fluxes, aparent magnitudees, and absolute magnitues of the stars.

An HR diagram is made of M(F336 - F555) by MF336 for units of both apparent and absolute magnitude. 

Histograms of the FWHM of detected peaks are also made. 

# Contents:
1. "data/" includes the folders F336W and F555W which include 3 fits files of images taken in theose respective filters (taken from brightspace).
2. optical_detectors.py: the script that will do the assignment
3. assignemnt_notebook.pynb: the jupyter notebook used to complete the assignment which shows where I tested parts of the code
4. script_testing.pynb: The notebook where I made the optical_detectors.py script and tested it using %run to make sure it runs with no string or plotting errors

# Notes:

When I run my script, I get a warning that says covarience could not be estimated at some point, but the code does then run 
(For context, this may be happening when I could not get a good gaussian fit for a peak, and in that case the code should just add it to a list of non-stars and move on). 
Hopefully that means that it will also run when anyone else uses it on their terminal. 
