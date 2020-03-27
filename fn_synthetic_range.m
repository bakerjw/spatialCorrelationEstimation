function [h, synthetic] = fn_synthetic_range (station_lat, station_long, recIdx, recsPerEQ, eventIdx, EQ_name_string, options)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Jack Baker 2/1/2019
% last updated 12/28/2019
%
% Simulate synthetic data with a given station configuration and assumed 
% semivariogram, and then re-estimate the semivariogram. 
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


% initialize matrices of simulated ranges and sills
for k=1:length(options.fitMethod)
    synthetic.sill{k} = ones(length(recIdx),options.numSims); 
    synthetic.range{k} = zeros(length(recIdx),options.numSims);
end

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
    %RHO = exp(-distance .^ 0.55 ./ options.syntheticRange); 
    
    % Simulate data and estimate variograms
    for j = 1:options.numSims
        rng(j);
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
        
        % perform variogram fits, using several techniques
        [sill, range] = fit_vario(h, gamma, nPairs, options);
        
        % copy fitted results into separate matrices for each fitting technique, to ease later analyses
        for k=1:length(options.fitMethod)
            synthetic.range{k}(i,j) = range(k);
            % note, if sills are being fitted, add code here to rearrange the results as well
        end
    
    end
    
    
    % Plot variogram replicates for this earthquake, if desired
    if options.plotFig
        hPlot = 0:0.5:options.maxR; % use a finer distance resolution for plotting semivariograms
        
        hf = figure;
        set(hf, 'Visible', 'off'); % don't show the figure on the screen
        h1=plot(h,1-exp(-3.*h./options.syntheticRange), '-k', 'linewidth', 2);
        hold on
        for j = 1:options.numSims
            h2=plot(hPlot,synthetic.sill{1}(i,j)*(1-exp(-3.*hPlot./synthetic.range{1}(i,j))), 'Color',[0.7 0.7 0.7], 'linewidth', 1);
        end
        legend([h1 h2],'Assumed model','Fitted replicates', 'location', 'southeast')
        xlabel('h (km)');
        ylabel('\gamma(h)');
        set(gca, 'xlim', [0 80])
        FormatFigureBook
        title(EQ_name_string{i})
        FormatFigureBook
        print('-dpdf', [options.figurePath 'perEvent/' num2str(i) '_semivariogramReplicates.pdf']); % save the figure to a file
        close all

    end
end




end