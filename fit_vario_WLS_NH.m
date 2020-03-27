function [sill, range, errors] = fit_vario_WLS_NH(h, gamma, nPairs, options)

% fit variogram sill and range using a weighted least-squares opimization criteria
% Created by Jack Baker 12/28/2019
% modified by Yilin Chen 02/11/2020
%
% The weight function is W = n/h.^2, per the default variogram fittting
% function in the R gstat package (https://github.com/r-spatial/gstat/)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input Variables
%   h = vector of separation distances (lags)
%   gamma = empirical variogram
%   nPairs = number of data pairs at each distance
%   options.fixedSill =1 to pre-assume a sill of 1, =0 to fit sill from
%           data -- not implemented in this function
%   options.funcForm = 1 for sill*(1-exp(-3 * h / range)); 
%                    = 2 for sill*(1-exp(-h^0.55 / range))
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

W = nPairs ./ (h.^2); % weighting of errors when optimizing to find the range
vario_fun = vario_fun_form(options);

if options.fixedSill
    sill = 1;
    
    % Exhaustive search option
    rangeVals = 1:0.2:120;    
    for i=1:length(rangeVals)
        score(i) = sum(W .*  (gamma -  vario_fun([rangeVals(i)], h) ).^2);
    end
    % figure
    % plot(rangeVals,score)
    
    [~,idx]= min(score);
    range=rangeVals(idx);
    errors = W / (sum(W) / length(W)) .*  (gamma -  vario_fun([range], h) ).^2; % compute error values for evaluation of the algorithm

else
    x0 = [1 25]; % starting value of sill and range
    fun = @(b,x) b(1)*(1 - exp(-3 .* x ./ b(2))); % function to be optimized
    nlm = fitnlm(h,gamma,fun,x0,'Weight', W); % perform weighted least-squares curve fit
    
    % rename variables to pass back out of function
    sill = nlm.Coefficients{1, 1};
    range = nlm.Coefficients{2, 1};
end

end





 
        
