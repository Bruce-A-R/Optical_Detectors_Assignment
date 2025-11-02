# Optical_Detectors_Assignment
This repository is my optical detectors assignment.

When the script optical_detectors.py is run, you can add arguments for the path to folders F336W and F555W (in that order). It should run fine if you do not, because
if paths are not specified the script should assume they are in the "data" folder of this github, and they are. I just personally think it makes more sense to specify your file path each time. 

# What the script will do

Using data from the F336W and F555W filter images, this script will combine  will create a catalog of stars with shared coordinates across both 

# Contents:
1. "data/" includes the folders F336W and F555W which include 3 fits files of images taken in theose respective filters (taken from brightspace).
2. optical_detectors.py: the script that will do the assignment
3. assignemnt_notebook.pynb: the jupyter notebook used to complete the assignment which shows where I tested parts of the code
4. script_testing.pynb: The notebook where I made the optical_detectors.py script and tested it using %run to make sure it runs with no string or plotting errors

# Notes:

When I run my script, I get a warning that says covarience could not be estimated at some point, but the code does then run 
(For context, this may be happening when I could not get a good gaussian fit for a peak, and in that case the code should just add it to a list of non-stars and move on). 
Hopefully that means that it will also run when anyone else uses it on their terminal. 

My catalog is also not printing out all the collumns in the terminal sadly, so I will include a csv version of it in the .zip file with my results incase that is necessary. 
