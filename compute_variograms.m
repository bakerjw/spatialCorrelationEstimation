% Main script to perform analysis and make figures 
% Jack Baker
% Last updated 10/2/2019 to incorporate new analysis
% Updated 12/19/2019 to modularize calculations
% Updated 3/17/2020 to better organize permutations of semivariogram fitting

clear; close all; clc;

%% user inputs

options.figurePath = 'figures/'; % location to print figures
if ~isfolder(options.figurePath)
    mkdir(options.figurePath)
end

options.TStar = 1; % spectral acceleration period of interest 
options.minNumRecordings = 40; % how many recordings are needed for an earthquake to be considered
options.fixedSill = 1; % =1 to pre-assume a sill of 1, =0 to fit sill from data (implemented only for some techniques)
options.renormalize = 1; % =1 to renormalize each event's data to have a standard deviation of 1, =0 otherwise
options.funcForm = 1; % =1 for exponential semivariogram, = 2 for exp(-a * h^0.55)
options.WLScoeff = 5; % coefficient for fit_vario_WLS weighting
options.eventNums = [7 30]; % two example events to plot histograms

% semivariogram parameters
options.maxR = 60; % maximum considered distance when computing semivariograms
options.binSize = 1; % 
options.plotFig = 0; % plot diagnostic figures for each earthquake
options.plotFigExtra = 0; % plot extra figures for each earthquake

% synthetic replicates parameters
options.numSims = 100;
options.syntheticRange = 30;
options.funcForm = 1;

% which fitting techniques to use
options.fitMethod = 1:6;

% labels for lines
options.linespec{1} = '--k';
options.linespec{2} = '--c';
options.linespec{3} = '-c';
options.linespec{4} = '-g';
options.linespec{5} = '-b';
options.linespec{6} = '-k';



%% load residuals data
addpath('data/')
load NGA_W2_corr_meta_data %Sa_RotD50
load NGA_W2_meta_data station_long station_lat hypo_long hypo_lat %EQ_name EQ_year magnitude closest_D
load mixedEffectsResids resid_RotD50Within medianPred

residData = resid_RotD50Within; % consider within-event residuals
tIdx = find(Periods == options.TStar); % get index for Period of interest



%% search for well-recorded events
for i=1:max(eqid)
    recsPerEQ{i} = find(eqid==i & ~isnan(residData(:, tIdx)));
    numRecs(i) = length(recsPerEQ{i});
    recsPerEQNoFilter{i} = find(eqid==i); % all recordings from earthquake, regardless of whether usable
    numRecsNoFilter(i) = length(recsPerEQNoFilter{i});
    if numRecs(i)>0
        exampleRec(i) = recsPerEQ{i}(1); % get the first record from this event (to use for finding EQ name, etc.)
    else
        exampleRec(i) = nan;
    end
end

% sort and trim event data
[sortedNumRecs,eventIdxAll] = sort(numRecs,'descend'); 
eventIdx = eventIdxAll(sortedNumRecs>=options.minNumRecordings)'; % select the events with a sufficient number of recordings, and transpose
recIdx = exampleRec(eventIdx);

eventIdxLg = eventIdxAll(sortedNumRecs>=100)'; % a second index of only very well recorded earthquakes
recIdxLg = exampleRec(eventIdxLg);


% format numerical-valued earthquake names to strings, for output to a table
for i=1:length(recIdx)
   temp = EQ_name{recIdx(i)};
   if isnumeric(temp)
       EQ_name_string{i,1} = num2str(temp);
   else
       EQ_name_string{i,1} = temp;
   end
end

% make a summary table of the events (formatted to be a Latex table)
summaryData = [EQ_name_string EQ_year(recIdx) num2cell(eventIdx) num2cell(magnitude(recIdx)) num2cell(numRecsNoFilter(eventIdx)')  num2cell(numRecs(eventIdx)') ]';
sprintf('%s & %.0f & %.0f & %.1f & %.0f & %.0f \\\\\n', summaryData{:})


%% collect data for each event, for more detailed analysis

for i=1:length(recIdx)
    idx{i} = recsPerEQ{eventIdx(i)}; % load allowable records for this event
    resids{i} = residData(idx{i}, tIdx);
    amplitudes{i} = Sa_RotD50(idx{i},tIdx);
    lats{i} = station_lat(idx{i});
    longs{i} = station_long(idx{i});
end


save('singleEventData.mat', 'options', 'resids', 'amplitudes', 'lats', 'longs')
 

% return


%% summary analysis for each event

nBins = round(options.maxR/options.binSize); % how many semivariogram bins will there be
gamma  = zeros(length(recIdx),nBins); % initialize matrix of semivariograms
nPairs = zeros(length(recIdx),nBins); % initialize matrix of number of station pairs


for i=1:length(recIdx)
    idx = recsPerEQ{eventIdx(i)}; % load allowable records for this event
    resids = residData(idx, tIdx);
    if options.renormalize
        resids = resids./std(resids); % renormalize residuals
    end
        
    % residuals vs distance
    b = regress(resids, [ones(size(closest_D(idx))) closest_D(idx)]);    % regression coefficients for distance vs resid
    distanceResidSlope(i) = b(2); % store regression slope for later analysis
    
    % semivariogram
    [sill(:,i), range(:,i), h, gamma(i,:), nPairs(i,:), methodName, methodNameShort] = fn_compute_variogram(station_lat(idx), station_long(idx), resids, options);
    %[sillL(i,1 ), rangeL(i,1), sillH(i,1), rangeH(i,1), sillLS(i,1), rangeLS(i,1), sillWLS(i,1), rangeWLS(i,1), h, gamma(i,:), nPairs(i,:)] = fn_compute_variogram (station_lat(idx), station_long(idx), resids, options);
    if options.plotFig % finish plot that was started inside fn_compute_variogram
        title(EQ_name_string{i})
        FormatFigureBook
        print('-dpdf', [options.figurePath 'perEvent/' num2str(i) '_exampleSemivariogram.pdf']); % save the figure to a file
        close all
    end

    if options.plotFigExtra
        % Plot residuals vs distance
        distRef = [0 300]; % reference distance values to use in residual plots
        
        hf = figure;
        set(hf, 'Visible', 'off'); % don't show the figure on the screen
        plot(closest_D(idx), resids, '.k')
        hold on
        plot(distRef, [0 0], '-k')
        plot(distRef, [b(1) + b(2).*distRef], '-b') % add a regression line
        xlabel('Closest distance (km)');
        ylabel('\delta W')
        title(EQ_name_string{i})
        legend('Data', '\delta W = 0', 'Regression line')
        set(gca, 'ylim', [-6 6])
        FormatFigureBook
        print('-dpdf', [options.figurePath 'perEvent/' num2str(i) '_residsVsDistance.pdf']); % save the figure to a file
        
        % Plot amplitudes versus distance
        hf = figure;
        set(hf, 'Visible', 'off'); % don't show the figure on the screen
        loglog(closest_D(idx), Sa_RotD50(idx,tIdx), '.')
        hold on        
        xlabel('Closest distance (km)');
        ylabel(['SA(' num2str(options.TStar) 's) [g]'])
        title(EQ_name_string{i})        
        FormatFigureBook
        print('-dpdf', [options.figurePath 'perEvent/' num2str(i) '_ampsVsDistance.pdf']); % save the figure to a file
         
        % Plot latitudes and longitudes
        hf = figure;
        set(hf, 'Visible', 'off'); % don't show the figure on the screen
        plot(station_long(idx), station_lat(idx), '.')
        hold on        
        plot(hypo_long(idx(1)), hypo_lat(idx(1)), 'pk')
        xlabel('Longitude');
        ylabel('Latitude')
        title(EQ_name_string{i})
        FormatFigureBook
        print('-dpdf', [options.figurePath 'perEvent/' num2str(i) '_LatLong.pdf']); % save the figure to a file      
        
        close all
    end
    
end

% save workspace to avoid later recomputation
save main_data station_lat station_long magnitude recIdx recIdxLg recsPerEQ eventIdx eventIdxLg EQ_year EQ_name_string options nPairs h gamma sill range methodName methodNameShort numRecs






