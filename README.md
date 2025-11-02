# Optical_Detectors_Assignment
This repository is my optical detectors assignment.

Contents:
1. "data" includes the folders F336W and F555W which include 3 fits files of images taken in theose respective filters (taken from brightspace).
2. optical_detectors.py: the script that will do the assignment
3. assignemnt_notebook.pynb: the jupyter notebook used to complete the assignment which shows where I tested parts of the code
4. script_testing.pynb: The notebook where I made the optical_detectors.py script and tested it using %run to make sure it runs with no string or plotting errors

Notes:

When I run my script, I get an error that says covarience could not be estimated at some point, but the code does then run 
(For context, this may ne happening when I could not get a good gaussian fit for a peak, and in that case the code should just add it to a list of non-stars and move on). 
Hopefully that means that it will also run when anyone else uses it on their terminal. 
