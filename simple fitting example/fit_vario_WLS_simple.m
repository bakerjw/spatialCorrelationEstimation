function [sill, range] = fit_vario_WLS_simple(h, gamma, nPairs, options)

% A simple function to fit a variogram sill and range using a weighted 
% least-squares opimization criteria
%
% Created by Yilin Chen
% last modified by Jack Baker, 3/27/2020
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input Variables
%   h = vector of separation distances (lags)
%   gamma = empirical variogram
%   nPairs = number of data pairs at each distance
%   options.WLScoeff = coefficient for weight taper
%
% Output Variables
%   sill = fitted sill for an exponential variogram model
%   range = fitted range for an exponential variogram model
%   errors = error value at each separation distance, as computed using this algorithm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% remove nan's from semivariogram
idx = find(~isnan(gamma));

gamma = gamma(idx);
h = h(idx);
nPairs = nPairs(idx);

% compute weights
W = nPairs .* exp(-h/options.WLScoeff); % use an exponentially decaying weighting

% get functional form
vario_fun = vario_fun_form();

x0 = [1 25]; % starting estimate of sill and range
nlm = fitnlm(h,gamma,vario_fun,x0,'Weight', W); % perform weighted least-squares curve fit

% rename variables to pass back out of function
sill = nlm.Coefficients{1, 1};
range = nlm.Coefficients{2, 1};

end





 function [fun] = vario_fun_form()
% variogram functional form to be fit
% The variogram form is exponential:  gamma = sill*(1-exp(-3*h/range))

fun = @(b,x)  b(1).*(1 - exp(-3 .* x ./ b(2)));
    
end
        
