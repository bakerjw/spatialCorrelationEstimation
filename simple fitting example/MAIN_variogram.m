% A simple script to compute a semivariogram and fit a model
%
% Created by Jack Baker, 3/27/2020


% user inputs
options.renormalize = 1; % =1 to renormalize each event's data to have a standard deviation of 1, =0 otherwise
options.WLScoeff = 5; % coefficient for fit_vario_WLS weighting
options.maxR = 60; % maximum considered distance when computing semivariograms
options.binSize = 1;  

% load residuals data
load example_data % load some example data to analyze. 
    % lats = vector of latitudes of the observations
    % longs = vector of longitudes of the observations
    % resids = vector of values to be analyze
    

% compute semivariogram
[sill, range, h, gamma] = fn_simple_variogram(lats, longs, resids, options);





