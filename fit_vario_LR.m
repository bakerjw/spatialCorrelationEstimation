function [sill, range, errors] = fit_vario_LR(h, gamma, options)

% fit variogram using the approach with a transformation to allow linear
% regression estimation
%
% For reference, see the following document:
% 
% Loth, C., and Baker, J. W. (2013). "A spatial cross-correlation model for 
% ground motion spectral accelerations at multiple periods." Earthquake 
% Engineering & Structural Dynamics, 42(3), 397-417.
%
% Created by Christophe Loth, modified by Jack Baker 10/2/2019
% heavily modified by Yilin Chen 02/11/2020
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input Variables
%   h = vector of separation distances (lags)
%   gamma = empirical variogram
%   options.fixedSill =1 to pre-assume a sill of 1, =0 to fit sill from
%   data 
%   options.funcForm = 1 for sill*(1-exp(-3 * h / range)); 
%                    = 2 for sill*(1-exp(-h^0.55 / range))
% Output Variables
%   sill = fitted sill for an exponential variogram model (should be one in 
%   the case of a direct semivariogram (single random variable)
%   range = fitted range for an exponential variogram model
%   errors = error value at each separation distance, as computed using this algorithm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% remove nan's from semivariogram
idx = find(~isnan(gamma));
gamma = gamma(idx);
h = h(idx);


if options.fixedSill
    % assume the sill
    sill = 1;
else
    % Estimate the sill
    y=0:0.005:1.5;
    kernel=zeros(length(y),1);
    sigma=0.1; % standard deviation of the kernel (user specified)
    
    for i=1:length(y)
        kernel(i)=sum(exp(-(gamma-y(i)).^2/sigma));
    end
    
    sill=y(kernel==max(kernel));
    
end

inv_h=1./(h+0.01); % weight the observations inverse to distance

vario_fun = vario_fun_form(options);

ydata = log(sill-min(sill-0.01,gamma));

rangeVals = 1:0.2:120;

for i=1:length(rangeVals)
    ypred = log(sill-sill*vario_fun(rangeVals(i), h)); % biased
    %ypred = log(sill-min(sill-0.01, sill*vario_fun(rangeVals(i), h))); % unbiased
    score(i) = sum(inv_h .^2 / (sum(inv_h .^2) / length(inv_h)) .* (ydata - ypred).^2);
end

[errors,idx] = min(score);
range = rangeVals(idx);

end





 
        
