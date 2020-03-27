function [sill, range, errors] = fit_vario_Cressie(h, gamma, nPairs, options)

% fit variogram sill and range using a weighted least-squares opimization 
% criterion (Equation 23) from  
%
% Cressie, N. (1985). "Fitting variogram models by weighted least squares." 
% Journal of the International Association for Mathematical Geology, 17(5), 
% 563-586.
%
% Created by Jack Baker 10/7/2019
% heavily modified by Yilin Chen 02/11/2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input Variables
%   h = vector of separation distances (lags)
%   gamma = empirical variogram
%   nPairs = number of data pairs at each distance
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
nPairs = nPairs(idx);
h = h(idx);

sill = 1; % assumed for now

vario_fun = vario_fun_form(options);

rangeVals = 1:0.2:100;

for i=1:length(rangeVals)
    score(i) = sum(nPairs .* (gamma./vario_fun([rangeVals(i)], h) - 1).^2);
end

[~,idx] = min(score);
range = rangeVals(idx);

errors = nPairs / (sum(nPairs) / length(nPairs)) .* (gamma./vario_fun([range], h) - 1).^2; % compute error values for evaluation of the algorithm

end





 
        
