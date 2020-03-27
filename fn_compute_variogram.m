function [sill, range, h, gamma, nPairs, methodName, methodNameShort] = fn_compute_variogram (lats, longs, values, options)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initial code from Christophe Loth, 09/12/12
% heavily modified by Jack Baker 2/1/2019
% last updated 3/17/2020
%
% Calculate empirical semivariograms from data, using several fitting
% techniques
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input Variables
%   lats = list of site lats
%   longs = list of site longs
%   values = list of the random variable values at each site
%   options.maxR = maximum distance to which the variogram is computed
%   options.binSize = distance interval accounted for by each computed variogram 
%       value
%   options.plotFig = 1 to plot semivariogram, =0 to not
%
% Output Variables
%   h = vector of separation distances (lags)
%   gamma = empirical variogram
%   sill = fitted sill for an exponential variogram model
%   range = fitted range for an exponential variogram model
%   nPairs = number of station pairs with a given separation distance
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


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


%% compute semivariogram

% initialize variables
h = zeros(nmax,1);
gamma = zeros(nmax,1);     

% compute empirical semivariance
for i=1:nmax
    [site1, site2] = find(distance_ratio == i);
    nPairs(i,1) = length(site1);
    h(i,1)=options.binSize/2+(i-1)*options.binSize;
    gamma(i,1) = (1/(2*(nPairs(i,1))))*sum((values(site1)-values(site2)).*(values(site1)-values(site2)));
end

%% perform fits, using several techniques

[sill, range, ~, methodName, methodNameShort] = fit_vario(h, gamma, nPairs, options);

% Old version of fitting functions
% [sillL, rangeL] = fit_vario_Loth(h, gamma, options); % Fit using the approach of Loth and Baker 
% [sillH, rangeH] = fit_vario_Heresi_exp(h, gamma, options); % Fit using weighted least squares (Heresi and Miranda)
% % [sillLS, rangeLS] = fit_vario_LS(h, gamma, options); % Fit an exponential model
% [sillLS, rangeLS] = fit_vario_Cressie(h, gamma, nPairs, options); % Fit using weighted least squares (Cressie)
% [sillWLS, rangeWLS] = fit_vario_WLS(h, gamma, nPairs, options); % Fit using weighted least squares (distance)
% [sillNH,rangeNH] = fit_vario_NH(h, gamma, nPairs, options); % Fit using n-h weighting



% Plot the results
if options.plotFig
    hPlot = 0:0.5:options.maxR; % use a finer distance resolution for plotting semivariograms
    legendText{1} = 'Empirical semivariogram';

    hf = figure;
    set(hf, 'Visible', 'off'); % don't show the figure on the screen
    hEmp=plot(h, gamma, '.k', 'linewidth', 2);
    hold on
    for i = 1:length(options.fitMethod)
        h(i) = plot(hPlot, sill(i) * (1-exp(-3.*hPlot./range(i))),options.linespec{i});
        legendText{i+1} = methodName{i};
    end
%         h2=plot(hPlot, sillL*  (1-exp(-3.*hPlot./rangeL)), 'Color',[0.7 0.7 0.7], 'linewidth', 2);
%         h3=plot(hPlot, sillH*  (1-exp(-3.*hPlot./rangeH)), 'Color','r', 'linewidth', 2);
%         h4=plot(hPlot, sillLS* (1-exp(-3.*hPlot./rangeLS)), 'Color','b', 'linewidth', 2);
%         h5=plot(hPlot, sillWLS*(1-exp(-3.*hPlot./rangeWLS)), 'Color','g', 'linewidth', 2);
    legend(legendText, 'location', 'southeast')
    xlabel('h [km]');
    ylabel('\gamma(h)');
    set(gca, 'xlim', [0 80])
    FormatFigureBook
end

end