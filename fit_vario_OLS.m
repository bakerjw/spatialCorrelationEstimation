function [sill, range, errors] = fit_vario_OLS(h, gamma, options)

% fit variogram using a least-squares opimization criteria
%
% Created by Yilin Chen, modified by Jack Baker 10/2/2019
% heavily modified by Yilin Chen 02/11/2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input Variables
%   h = vector of separation distances (lags)
%   gamma = empirical variogram
%   options.fixedSill =1 to pre-assume a sill of 1, =0 to fit sill from data 
%   options.funcForm = 1 for sill*(1-exp(-3 * h / range)); 
%                    = 2 for sill*(1-exp(-h^0.55 / range))
% Output Variables
%   sill = fitted sill for an exponential variogram model
%   range = fitted range for an exponential variogram model
%   errors = error value at each separation distance, as computed using this algorithm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% remove nan's from semivariogram
idx = find(~isnan(gamma));
gamma = gamma(idx);
h = h(idx);

if options.fixedSill
    sill = 1;
    
    vario_fun = vario_fun_form(options);

    rangeVals = 1:0.2:120;    
    for i=1:length(rangeVals)
        score(i) = sum((gamma -  vario_fun([rangeVals(i)], h) ).^2);
    end
    
    [errors,idx]= min(score);
    range=rangeVals(idx);

else
    x0 = [1 25]; % starting value of sill and range
    lb = [0.1 2]; % lower bound value of sill and range
    ub = [4 200]; % upper bound value of sill and range (set large so that it shouldn't be a constraint)
    
    fun = @(x,h)(x(1)*(1 - exp(-3 .* h ./ x(2)))); % function to be optimized
    
    x = lsqcurvefit(fun, x0, h, gamma, lb, ub, 'Display', 'off'); % perform least-squares curve fit
    
    % rename variables to pass back out of function
    sill = x(1);
    range = x(2);
end

end





 
        
