function [sill, range, loss, methodName, methodNameShort] = fit_MLL(lats, longs, values, options)

% fit variogram sill and range using the log-likelihood of a multivariate 
% normal distribution as in GP regression.
% Created by Lukas Bodenmann 04/13/2022
% This method finds the variance and the lengthscale of the covariance 
% function that maximize the log-likelihood of the observations.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input Variables
%   lats = list of site lats
%   longs = list of site longs
%   values = list of the random variable values at each site
%   options.fixedSill   = 1 to pre-assume a sill of 1, 
%                       = 0 to fit sill from data
%   options.funcForm = 1 for sill*(1-exp(-3 * h / range)); 
%                    = 2 for sill*(1-exp(-h^0.55 / range))
%
% Output Variables
%   sill = fitted sill for an exponential variogram model
%   range = fitted range for an exponential variogram model
%   loss = Average log-likelihood of the fitted model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Assemble the input matrix
X = [lats longs];

% Treat records from stations that were closer than 5m as duplicates and 
% delete them
% This affects less than 10 station pairs over all events with more than 
% 40 recordigs.
c = 1;
idx_dupl = int16.empty();
dist = pos2distv(X, X);
for i =1:length(dist)
    p = sum(dist(i,:)<0.005);
    if p > 1
        idx_dupl(c) = i;
        c = c+1;
    end
end
lats(idx_dupl) = []; longs(idx_dupl) = [];
values(idx_dupl) = [];
X = [lats longs];

% Specify the covariance functions depending on the variogram forms
if options.funcForm == 1 && options.fixedSill == 1
    kfcn = @(XN, XM, theta) exp(-3 * pos2distv(XN,XM)/exp(theta(1)));
    theta0 = log([45]); % set initial value in log scale
    gprMdl = fitgp(X, values, kfcn, theta0);
    range = exp(gprMdl.KernelInformation.KernelParameters); % output
    sill = 1.0;
elseif options.funcForm == 2 && options.fixedSill == 1
    kfcn = @(XN, XM, theta) exp(-1 * pos2distv(XN,XM).^0.55 / exp(theta(1)));
    theta0 = log([45]);
    gprMdl = fitgp(X, values, kfcn, theta0);
    range = exp(gprMdl.KernelInformation.KernelParameters);
    sill = 1.0;
elseif options.funcForm == 1 && options.fixedSill == 0
    kfcn = @(XN, XM, theta) exp(theta(2)) * exp(-3 * pos2distv(XN,XM)/exp(theta(1)));
    theta0 = log([45, 1.0]);
    gprMdl = fitgp(X, values, kfcn, theta0);
    params = exp(gprMdl.KernelInformation.KernelParameters);
    range = params(1,1);
    sill = params(2,1);
elseif options.funcForm == 2 && options.fixedSill == 0
    kfcn = @(XN, XM, theta) exp(theta(2)) * exp(-1 * pos2distv(XN,XM).^0.55 / exp(theta(1)));
    theta0 = log([45, 1.0]);
    gprMdl = fitgp(X, values, kfcn, theta0);
    params = exp(gprMdl.KernelInformation.KernelParameters);
    range = params(1,1);
    sill = params(2,1);
end

loss = gprMdl.LogLikelihood/length(lats); % Normalized joint log-likelihood

methodName = 'GP Regression, Max Likelihood';
methodNameShort = 'GPR, ML';

% Wrapper for the matlab gp fit method
% A small noise value ('Sigma') is assigned for numerical stability
function GPmodel = fitgp(X, y, kernel, kernel_init)
GPmodel = fitrgp(X,y,'BasisFunction', 'none', ...
        'KernelFunction',kernel,'KernelParameters',kernel_init, ...
        'FitMethod','exact','Optimizer','lbfgs', ...
        'Sigma',1e-4, 'ConstantSigma', true, 'SigmaLowerBound', 1e-6);
end

end