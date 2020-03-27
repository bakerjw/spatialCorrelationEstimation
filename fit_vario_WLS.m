function [sill, range, errors] = fit_vario_WLS(h, gamma, nPairs, options)

% fit variogram sill and range using a weighted least-squares opimization criteria
% Created by Yilin Chen, modified by Jack Baker 12/28/2019
% heavily modified by Yilin Chen 02/11/2020
% The variogram form is exponential:  gamma = sill*(1-exp(-3*h/range))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input Variables
%   h = vector of separation distances (lags)
%   gamma = empirical variogram
%   nPairs = number of data pairs at each distance
%   options.fixedSill =1 to pre-assume a sill of 1, =0 to fit sill from
%   data
%   options.WLScoeff = coefficient for weight taper
%   options.funcForm = 1 for sill*(1-exp(-3 * h / range)); 
%                    = 2 for sill*(1-exp(-h^0.55 / range))
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


% W = exp(-options.WLScoeff * h); % use an exponentially decaying weighting, when optimizing to find the range

% weight option with number-of-pairs weighting
nPairs = nPairs(idx);
W = nPairs .* exp(-h/options.WLScoeff); % use an exponentially decaying weighting, when optimizing to find the range

% weight by number of stations and variogram value, https://pdf.sciencedirectassets.com/271720/1-s2.0-S0098300408X00068/1-s2.0-S0098300408000666/main.pdf?X-Amz-Security-Token=IQoJb3JpZ2luX2VjEEsaCXVzLWVhc3QtMSJFMEMCIFtpS1suV%2FMqW%2F5m5vn6uZVZGUyBKQNVoHheLBA2letsAh9%2Fk8oti9NcXjVmYeuwglbOudSX7q5cUHkEEdzj%2BKNUKr0DCLP%2F%2F%2F%2F%2F%2F%2F%2F%2F%2FwEQAhoMMDU5MDAzNTQ2ODY1Igyhj%2Bdujhn2Y2JABysqkQOzVAROp%2F3Jng7ltqHHclq7uWvZyNqIPOiduPLg9tJYVhDgD6I2nEjmnbqn%2F0%2BG9hTaYLoIcXFXS0yPSVyltlH3q9Li7IGC6o1VMJLHTGh22az73kKHCI4AXujypHDCm4vHSLYqisFI58U7TUgkISWkB6VQV7XgASFC1XUAgcg5LbjEu1e7t8N%2Fo4lVpxsnQalhv%2FHPOBZp0xxqFoA1YvK4LJpnD6lkEdG2bxkHLaGt4d8m5ZFWPQi9uwaeuCC5Z1lJEEYO%2FRjrdXh8qaEDq%2FWZH%2BaFIbN%2B%2ByO89eYRnEfKHtK2NkdBaE9CIrNqRV8Wg64rLxGKfpvZSJsQEwH0nqua3dKzVFQPsvVwzw9hUSRCtGnsmLAHE749VgEtXYdXW4DMbhHc81YXK6Bjea1slmrtMBPyvTY8umfGCJE9olFVF2AmWqdhLzEH6D7OGdxK0Mqa1IoK3Ncv7kHqF37ybUzPFq4qn323MyXSvXV7eapvfdOibPXt1ZjRUCACUIumueE83m3KeJysvRGaxGoo1LqWcTDQmqDwBTrtAexif2dSyGs52VvaM%2BjYCycCQqDKWZskkKQ7aeiv0OcMKnxP7b39FzQKCq6ZrPejJUccXdHz2fcjRNcZYnR3sGb%2BuU3rEIZCK8eoqebgvzQ9bfgPnHfLaGjjP9B1YVM0uH99o4wR01mIPdewpTfLpZnadMuw%2BPlFSEl%2FBpiTAfXwz8X%2FVN5Qv0lzYCw18mqx8B%2FfO%2Blcmg33GgybXkongukrY6p2770irDlF8PVPqgRSXBR0PAUkvDiDKEm3qQOodomGnoUXFRZ%2F7aopYQ%2BMNfEM4QZSMtcAqQ7pZUjeJXD27ou%2B7kxFxHjNxfC5Tg%3D%3D&X-Amz-Algorithm=AWS4-HMAC-SHA256&X-Amz-Date=20191229T034142Z&X-Amz-SignedHeaders=host&X-Amz-Expires=300&X-Amz-Credential=ASIAQ3PHCVTYQBWHVYKM%2F20191229%2Fus-east-1%2Fs3%2Faws4_request&X-Amz-Signature=7168d19608cedfae29f0a2797306ca331daf18a9bbcf7539542d5e24b40ce36c&hash=950464174b7f78c9b570090cbb056f7770e4ad5660260afa47252bff90a33082&host=68042c943591013ac2b2430a89b270f6af2c76d8dfd086a07176afe7c76c2c61&pii=S0098300408000666&tid=spdf-45f699a2-f34d-497c-b9d5-5ca6578cd6c1&sid=80dd761f5d50864aba78cab3fb82d8836d68gxrqa&type=client
% W = nPairs ./ (gamma.^2); 
% W = nPairs ./ (gamma); 
% W = nPairs ./ (gamma); 

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
    nlm = fitnlm(h,gamma,vario_fun,x0,'Weight', W); % perform weighted least-squares curve fit
    % rename variables to pass back out of function
    sill = nlm.Coefficients{1, 1};
    range = nlm.Coefficients{2, 1}; 
    errors = W .*  (gamma - vario_fun([sill range], h) ).^2; % compute error values for evaluation of the algorithm
end
end





 
        
