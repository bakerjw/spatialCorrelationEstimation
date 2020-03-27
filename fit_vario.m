function [sill, range, errors, methodName, methodNameShort] = fit_vario(h, gamma, nPairs, options)
% fit variogram sill and range using various methods  
%
% Created by Yilin Chen 02/22/2020
% Modified by Jack Baker, 3/17/2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input Variables
%   h = vector of separation distances (lags)
%   gamma = empirical variogram
%   nPairs = number of data pairs at each distance
%   options.fitMethod = vector of fitting methods of interest
%                     = 1 Ordinary Least Squares
%                     = 2 Weighted Least Squares with n*exp(-distance) weighting
%                     = 3 Weighted Least Squares with n/h^2 weighting 
%                     = 4 Cressie (1985)
%                     = 5 Fisher Transformation (Heresi and Miranda 2019) 
%                     = 6 Linear Regression transformation (see Loth and Baker 2013)
% Output Variables
%   sill = fitted sill for an exponential variogram model
%   range = fitted range for an exponential variogram model
%   errors = error value at each separation distance, as computed using
%   different algorithms
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = 1:length(options.fitMethod) % for each fitting method of interest
    method = options.fitMethod(i);
    
    switch method
        case 1
            [sill(i), range(i), errors{i}] = fit_vario_OLS(h, gamma, options);
            methodName{i} = 'Ordinary Least Squares';
            methodNameShort{i} = 'OLS';
            
        case 2
            [sill(i), range(i), errors{i}] = fit_vario_WLS(h, gamma, nPairs, options);
            methodName{i} = 'Weighted Least Squares';
            methodNameShort{i} = 'WLS';
            
        case 3
            [sill(i), range(i), errors{i}] = fit_vario_WLS_NH(h, gamma, nPairs, options);
            methodName{i} = 'Weighted Least Squares, w = n/h^2';
            methodNameShort{i} = 'WLS, nh';

        case 4
            [sill(i), range(i), errors{i}] = fit_vario_Cressie(h, gamma, nPairs, options);
            methodName{i} = 'Cressie';
            methodNameShort{i} = 'CR';

        case 5
            [sill(i), range(i), errors{i}] = fit_vario_FT(h, gamma, options);
            methodName{i} = 'Fisher Transform';
            methodNameShort{i} = 'FT';

        case 6
            [sill(i), range(i), errors{i}] = fit_vario_LR(h, gamma, options);
            methodName{i} = 'Linear Regression';
            methodNameShort{i} = 'LR';
            
            
            
    end
    
end    
    

    
end