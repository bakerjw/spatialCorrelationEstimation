# Spatial correlation estimation

This repository provides Matlab scripts to model spatial correlation in earthquake ground motions, and quantify estimation uncertainty. Further documentation can be found in the following paper:

Baker, Jack. W., and Chen, Yilin. (2020). “Ground motion spatial correlation fitting methods and estimation uncertainty.” Earthquake Engineering & Structural Dynamics, 49(15), 1662–1681. 
https://www.jackwbaker.com/Publications/Baker_Chen_(2020)_Spatial_fitting,_EESD.pdf, 
https://doi.org/10.1002/eqe.3322

The script MAIN_script.m in the root directory will produce the data and figures used in the above manuscript.

The script MAIN_variogram.m in the "simple fitting example" folder will perform a very simple computation of a semivariogram and model fit, to illustrate the basic mechanics of geostatical calculations.

The NGA-West2 data were updated in March 2022 to correct errors in the Z1 and region values associated with the ground motion database, and the code now produces figures that are slightly updated relative to the originally published results. Thank you to Eduardo Miranda and Alan Poulos for identifying the errors and helping resolve them.
