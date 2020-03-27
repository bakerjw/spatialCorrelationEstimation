function [synthetic, WLSnh] = fn_WLS_test (station_lat, station_long, recIdx, recsPerEQ, eventIdx, options, WLScoeffVals)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Jack Baker 2/1/2019
% last updated 3/17/2020 to include a test of the n/h^2 WLS approach
%
% Simulate synthetic data with a given station configuration and assumed 
% semivariogram, to test performance of WLS function 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input Variables
%   station_lat = list of site lats
%   station_long = list of site longs
%   recIdx
%   options.numSims = number of replicates desired
%   options.maxR = maximum distance to which the variogram is computed
%   options.binSize = distance interval accounted for by each computed variogram 
%       value
%   options.syntheticRange = assumed range used to simulate synthetic data
%   options.plotFig = 1 to plot semivariogram, =0 to not
%   options.figurePath = path to folder to save figures
%
% Output Variables
%   h = vector of separation distances (lags)
%   synthetic = structure containting replicates of sills and ranges
%   estimated using various methods
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



for i=1:length(recIdx)
    idx = recsPerEQ{eventIdx(i)}; % load allowable records for this event
    
    % locations of considered stations for this event
    lats = station_lat(idx);
    longs = station_long(idx);
        
    % Compute site-site distances
    nsites=length(lats);
    nmax = round(options.maxR/options.binSize);
    distance = zeros(nsites);
    distance_ratio = zeros(nsites);
    for k=1:nsites
        for j=1:nsites
            distance(k,j) = pos2dist(lats(k),longs(k),lats(j),longs(j),1);
            distance_ratio(k,j) = round(distance (k,j)/options.binSize); % To group stations in bins
        end
    end
    % Compute correlation matrix, using assumed range
    RHO = exp(-3.*distance./options.syntheticRange); 
    
    % Simulate data and estimate variograms
    for j = 1:options.numSims
        values = mvnrnd(zeros(1,nsites), RHO); % simulate residuals
        
        % compute semivariance from simulated residuals
        h = zeros(nmax,1); % initialize vector
        gamma = zeros(nmax,1); % initialize vector
        for k=1:nmax
            [site1, site2] = find(distance_ratio == k);
            nPairs(k,1) = length(site1);
            h(k,1)=options.binSize/2+(k-1)*options.binSize;
            gamma(k,1) = (1/(2*(nPairs(k,1))))*sum((values(site1)-values(site2)).*(values(site1)-values(site2)));
        end
        
        for g = 1:length(WLScoeffVals) % estimate variograms using each weighting function
            options.WLScoeff = WLScoeffVals(g);
            [synthetic.sill{g}(i,j), synthetic.range{g}(i,j)] = fit_vario_WLS(h, gamma, nPairs, options); 
        end
        
        % get estimates from n/h^2 WLS approach
        [WLSnh.sill(i,j), WLSnh.range(i,j)] = fit_vario_WLS_NH(h, gamma, nPairs, options); 
       
    end
    
end




end