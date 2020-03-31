% Calculations to decompose the semivariogram variance 
%
% Created by Jack Baker
% Last modified 3/24/2019


clear; close all; clc;

load synthetic_data h synthetic
load main_data nPairs numRecs recIdx methodName options eventIdx range recsPerEQ



% focus on Weighted Least Squares results below
k = 2; % WLS index


%% statistics for each fitting method and earthquake

for i=1:length(recIdx)
    % compute mean range for each event
    meanMatrix(i) = mean(synthetic.range{k}(i,:));
    
    % compute standard deviation for each event
    stdMatrix(i) = std(synthetic.range{k}(i,:));
end



%% get estimation standard deviations and biases for well-recorded and poorly recorded earthquakes

nLg = 130; % how many recordings to fall in the 'large' category
nSm = 65; % how many recordings to fall in the 'small' category

idxLg = find(numRecs(eventIdx) > nLg);
estSigmaLg = mean(stdMatrix(idxLg))'; % note, we can take the average of the standard deviations rather than a pooled standard deviation, since there are the same number of data points for each earthquake

idxMd = find((numRecs(eventIdx) <= nLg) & (numRecs(eventIdx) > nSm) );
estSigmaMd = mean(stdMatrix(idxMd))'; 

idxSm = find(numRecs(eventIdx) <= nSm);
estSigmaSm = mean(stdMatrix(idxSm))'; 



%% get real data standard deviations

totalSigmaLg = std(range(k, idxLg)')';
totalSigmaMd = std(range(k, idxMd)')';
totalSigmaSm = std(range(k, idxSm)')'; 


sigmaRealLg = sqrt(totalSigmaLg.^2 - estSigmaLg.^2);
sigmaRealMd = sqrt(totalSigmaMd.^2 - estSigmaMd.^2);
sigmaRealSm = sqrt(totalSigmaSm.^2 - estSigmaSm.^2);



sprintf('$n_{EQ}$    %.1f      %.1f      %.1f  \\\\', [length(idxSm); length(idxMd); length(idxLg)])
sprintf('Total sigma    %.1f      %.1f      %.1f  \\\\', [totalSigmaSm; totalSigmaMd; totalSigmaLg])
sprintf('Estimation sigma    %.1f      %.1f      %.1f  \\\\', [estSigmaSm; estSigmaMd; estSigmaLg])
sprintf('Real sigma    %.1f      %.1f      %.1f  \\\\', [sigmaRealSm; sigmaRealMd; sigmaRealLg])

%% repeat calcs at multiple periods

if 1==1 % run new analysis
    options.fitMethod = k; % consider the same fitting method as specified above
    
    % load residuals data
    addpath('data/')
    load NGA_W2_corr_meta_data Periods eqid
    load NGA_W2_meta_data station_long station_lat  %EQ_name EQ_year magnitude closest_D
    load mixedEffectsResids resid_RotD50Within
    residData = resid_RotD50Within; % consider within-event residuals
    Periods(75:end) = []; % omit long period data with more limited results
    
    
    clear range sill % clear old variables because we are computing new range data for all periods
    
    for j = 1:length(Periods)
        options.TStar = Periods(j); % update this just in case it influences a later calculation
        
        for i=1:length(recIdx) % for each earthquaake
            recsPerEQ{i} = find(eqid==i & ~isnan(residData(:, j)));
            
            idx = recsPerEQ{eventIdx(i)}; % load allowable records for this event
            resids = residData(idx, j);
            if options.renormalize
                resids = resids./std(resids); % renormalize residuals
            end
            
            % semivariogram
            [~, range(j,i)] = fn_compute_variogram(station_lat(idx), station_long(idx), resids, options);
        end
    end
    

    save allPeriodsRanges range Periods

else % load previously calculated data
    load allPeriodsRanges
    
end


for j = 1:length(Periods)
    totalSigmaLg(j) = std(range(j, idxLg));
    totalSigmaMd(j) = std(range(j, idxMd));
    totalSigmaSm(j) = std(range(j, idxSm));
end
sigmaRealLg = sqrt(totalSigmaLg.^2 - estSigmaLg.^2);
sigmaRealMd = sqrt(totalSigmaMd.^2 - estSigmaMd.^2);
sigmaRealSm = sqrt(totalSigmaSm.^2 - estSigmaSm.^2);


%% plots

nAvg = 9; % optional parameter to smooth the results for ease of reading. nAvg=1 means no smoothing

figure
semilogx(Periods, movmean(totalSigmaLg,nAvg), '-k')
hold on
plot(Periods, movmean(totalSigmaMd,nAvg), '-b')
plot(Periods, movmean(totalSigmaSm,nAvg), '-r')

% semilogx(Periods, estSigmaLg*ones(size(Periods)), '--k')
% hold on
% plot(Periods, estSigmaMd*ones(size(Periods)), '--b')
% plot(Periods, estSigmaSm*ones(size(Periods)), '--r')

h1 = semilogx(Periods, movmean(sigmaRealLg,nAvg), '-k', 'linewidth', 2);
hold on
h2 = semilogx(Periods, movmean(sigmaRealMd,nAvg), '-b', 'linewidth', 2);
h3 = plot(Periods, movmean(sigmaRealSm,nAvg), '-r', 'linewidth', 2);
legend([h1 h2 h3], ['Earthquakes with n \geq ' num2str(nLg) ' recordings'], ['Earthquakes with ' num2str(nSm) ' \leq n < ' num2str(nLg) ' recordings'], ['Earthquakes with n < ' num2str(nSm) ' recordings'], 'location', 'Northeast')
set(gca, 'ylim', [0 50])
set(gca, 'xtick', [0.01 0.1 1])
set(gca, 'xticklabel', [0.01 0.1 1.0])
xlabel('Period [s]')
ylabel('Standard deviatiation')
FormatFigureBook
print('-dpdf', [options.figurePath 'variance_decomposition.pdf']); % save the figure to a file

