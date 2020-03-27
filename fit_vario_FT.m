function [sill, range, errors] = fit_vario_FT(h, gamma, options)

% fit variogram using the approach with a Fisher Transform
%
% For reference, see the following document:
% 
% Heresi, P., and Miranda, E. (2019). "Uncertainty in intraevent spatial 
% correlation of elastic pseudo-acceleration spectral ordinates." 
% Bulletin of Earthquake Engineering, 17(3), 1099-1115.
%
% Created by Jack Baker 10/8/2019
% heavily modified by Yilin Chen 02/11/2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input Variables
%   h = vector of separation distances (lags)
%   gamma = empirical variogram
%   options.fixedSill =1 to pre-assume a sill of 1, =0 to fit sill from
%           data -- not implemented in this function
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

sill = 1; % assume a sill of 1 (in terms of Heresi formulation, assume residuals have sigma=1)

rho = 1 - gamma;
ydata = 0.5*log( (1+rho)./(1-rho)); % perform a Fisher transformation on the data

vario_fun = vario_fun_form(options);

rangeVals = 1:0.2:120;

for i=1:length(rangeVals)
    gammapred = vario_fun([rangeVals(i)], h);
    ypred = 0.5*log( (1 + (1 - gammapred)) ./ (1 - (1 - gammapred)));
    score(i) = sum((ydata - ypred).^2);
end

[errors,idx] = min(score);
range = rangeVals(idx);

end






 
        
