function [sill, range, h, gamma, nPairs] = fn_simple_variogram (lats, longs, values, options)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initial code from Christophe Loth, 09/12/12
% last updated by Jack Baker 3/27/2020
%
% Calculate an empirical semivariograms from data, and fit a model
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input Variables
%   lats = list of site lats
%   longs = list of site longs
%   values = list of the random variable values at each site
%   options.maxR = maximum distance to which the variogram is computed
%   options.binSize = distance interval accounted for by each computed variogram 
%       value
%   options.renormalize = 1 to renormalize values and get a sill of 1
%
% Output Variables
%   sill = fitted sill for an exponential variogram model
%   range = fitted range for an exponential variogram model
%   h = vector of separation distances (lags)
%   gamma = empirical variogram
%   nPairs = number of station pairs with a given separation distance
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if options.renormalize
    values = values./std(values); % renormalize residuals, if desired
end

% Compute site-to-site distances
nsites=length(lats);
nmax = round(options.maxR/options.binSize);
        
distance = zeros(nsites);
distance_ratio = zeros(nsites);
for i=1:nsites
    for j=1:nsites
        distance(i,j) = pos2dist(lats(i),longs(i),lats(j),longs(j),1); 
        distance_ratio(i,j) = round(distance(i,j)/options.binSize); % To group stations in bins
    end
end

% compute empirical semivariance
h = zeros(nmax,1); % initialize variable
gamma = zeros(nmax,1); % initialize variables
nPairs = zeros(nmax,1); % initialize variables
for i=1:nmax
    [site1, site2] = find(distance_ratio == i);
    nPairs(i,1) = length(site1);
    h(i,1)=options.binSize/2+(i-1)*options.binSize;
    gamma(i,1) = (1/(2*(nPairs(i,1))))*sum((values(site1)-values(site2)).*(values(site1)-values(site2)));
end

% fit model
[sill, range] = fit_vario_WLS_simple(h, gamma, nPairs, options);

% Plot the results
hPlot = 0:0.5:options.maxR; % use a finer distance resolution for plotting semivariograms

figure;
plot(h, gamma, '.k', 'linewidth', 2);
hold on
plot(hPlot, sill * (1-exp(-3.*hPlot./range)));
legend('Empirical semivariogram', 'Fitted model', 'location', 'southeast')
xlabel('h [km]');
ylabel('\gamma(h)');
set(gca, 'xlim', [0 options.maxR])

end