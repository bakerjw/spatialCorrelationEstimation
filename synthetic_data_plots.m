% plots of results from synthetic data, to diagnose performance of the
% algorithms
%
% Created by Jack Baker, 10/3/2019
% Modified 3/17/2020 to use new fitting method format


clear; close all; clc;

load synthetic_data h synthetic
load main_data nPairs numRecs magnitude recIdx methodName options eventIdx range



figLabel = {'(a)', '(b)', '(c)', '(d)', '(e)', '(f)'};

yL = [0 120]; % y axis limits for some plots (temporarily hard coded)

%% statistics for each fitting method and earthquake

for i=1:length(recIdx)
    for k=1:length(options.fitMethod)
        % compute mean range for each event and fitting approach
        meanMatrix(i,k) = mean(synthetic.range{k}(i,:));

        % compute standard deviation for each event and fitting approach
        stdMatrix(i,k) = std(synthetic.range{k}(i,:));
        
        % compute mean square error for each event and fitting approach
        mseMatrix(i,k) = (meanMatrix(i,k)-options.syntheticRange).^2 + (stdMatrix(i,k)).^2;
    end
end

%% tabulate some standard deviations
nLg = 130; % how many recordings to fall in the 'large' category
nSm = 65; % how many recordings to fall in the 'small' category

idx{1} = find(numRecs(eventIdx) <= nSm);
idx{2} = find((numRecs(eventIdx) <= nLg) & (numRecs(eventIdx) > nSm) );
idx{3} = find(numRecs(eventIdx) > nLg);

idxVLg = find(numRecs(eventIdx) > 200); % indices of very well recorded earthquakes (where bias stabilizes)

for k=1:6 % for each fitting method
    estBias(k,1) = mean(meanMatrix(idxVLg,k)) - options.syntheticRange; % bias for 'large' category
    
    for i=1:3 % for each data subset
        estSigma(k,i) = mean(stdMatrix(idx{i},k))'; % note, we can take the average of the standard deviations rather than a pooled standard deviation, since there are the same number of data points for each earthquake
    end
end




summaryData = [methodName' num2cell(estBias) num2cell(estSigma) ]';
sprintf('%s  & %.1f  & %.1f & %.1f & %.1f \\\\\n', summaryData{:})

% Header line: \textbf{Event name} & \textbf{Year} & \textbf{Event index} & \textbf{Magnitude}  & \textbf{\# of usable stations} & \textbf{Estimation $\sigma$} & \textbf{Posterior mean $r$ [km]}         \\ 



%% plot two example histograms of ranges for illustration
k = 2; % pick a single fitting method to illustrate results

edges = [0:5:80];
yMax = 35; % vertical limit of plot

exampleStds = std(synthetic.range{k}(options.eventNums,:)')
exampleMeans = mean(synthetic.range{k}(options.eventNums,:)')

figure
for j=1:2
    subplot(1,2,j)
    histogram(synthetic.range{k}(options.eventNums(j),:),edges)     
    hold on
    plot(options.syntheticRange*[1 1], [0 yMax], '-k', 'linewidth', 2)
    set(gca, 'ylim', [0 yMax])
    xlabel('Range [km]')
    ylabel('Number of observations')
    legend('Histogram', 'Actual range')
    text(-0.1,-0.1,figLabel{j},'Units', 'Normalized', 'VerticalAlignment', 'Top')
end
Format2x1SubplotFigure
set(gcf, 'PaperSize', [7 2.5]); % make figure less tall than standard
set(gcf, 'PaperPosition', [0 0 7 2.5]);
print('-dpdf', [options.figurePath 'exampleHistograms.pdf']); % save the figure to a file



%% plots of ranges for each fitting method and replicate, versus number of stations
xL = [0 100*ceil(max(numRecs(eventIdx))/100)]; % x axis limits (go to the next round number above max value)
figure
for k=1:6
    subplot(2,3,k)
    hold on
    h2 = plot(repmat(numRecs(eventIdx)', 1, options.numSims/4), synthetic.range{k}(:,1:options.numSims/4), '.g'); % only plotting some of the simulations, to avoid clutter in the figure
    [xSort,meanSort] = fn_kernel_smooth(numRecs(eventIdx), meanMatrix(:,k));
    h3 = plot(xSort,meanSort, '-k', 'linewidth', 1);
    [xSort,stdSort] = fn_kernel_smooth(numRecs(eventIdx), stdMatrix(:,k));
    h4 = plot(xSort,meanSort+stdSort, '-.k', 'linewidth', 1);
    plot(xSort,meanSort-stdSort, '-.k', 'linewidth', 1)
    h1 = plot(xL, [1 1]*options.syntheticRange, ':k');

    if k==3 % add a legend on the second subplot
        legend([h1 h2(1) h3 h4], 'Correct range', 'Estimates', 'Moving Avg.', 'Moving Avg. \pm \sigma', 'location', 'Northeast')
    end
    
    axis([xL yL])
    xlabel(['Number of stations']);
    ylabel('Simulated ranges')
    title(methodName{k})
    text(-0.1,-0.1,figLabel{k},'Units', 'Normalized', 'VerticalAlignment', 'Top')
end
Format2x2SubplotFigure
print('-dpdf', [options.figurePath 'rangeReplicatesVsNumStations.pdf']); % save the figure to a file








