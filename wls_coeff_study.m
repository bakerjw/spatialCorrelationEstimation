% Use coordinates from real earthquake stations but generate synthetic
% residuals, to analyze performance of weighted least squares fitting 
% with various parameter values
%
% Created by Jack Baker, 12/28/2019
% Edited 3/17/2020

clear; close all; clc;

load main_data station_lat station_long recIdx recsPerEQ eventIdx options numRecs
figLabel = {'(a)', '(b)', '(c)', '(d)', '(e)', '(f)'};

%% which earthquakes to evaluate

eventIdxTest = eventIdx; % use all considered data
recIdxTest = recIdx; % use all considered data


%% collect estimated ranges for a set of Weighted Least Squares fitting approaches 

% WLS coefficient values to consider

if 1==1 % perform simulations
    WLScoeffVals = [0.5:0.5:6 7:15];
    rangeVals = [40:-5:20];
    
    for k = 1:length(rangeVals)
        options.syntheticRange = rangeVals(k);
        [WLSsynthetic{k}, WLSnh{k}] = fn_WLS_test(station_lat, station_long, recIdxTest, recsPerEQ, eventIdxTest, options, WLScoeffVals);
    end
   
    save WLS_coeff_data_multi_range WLSsynthetic WLScoeffVals WLSnh rangeVals
else % load previous data
    load WLS_coeff_data_multi_range
end

    
%% summary statistics

for k = 1:length(rangeVals)
    for j=1:length(WLScoeffVals)
        for i=1:length(recIdxTest)
            meanMatrix{k}(i,j) = mean(WLSsynthetic{k}.range{j}(i,:)); % mean range for each event and fitting approach
            stdMatrix{k}(i,j) = std(WLSsynthetic{k}.range{j}(i,:)); % standard deviation for each event and fitting approach
        end
        BIAS{k} = mean (abs(meanMatrix{k}-rangeVals(k)));
        AvgCOV{k} = mean (stdMatrix{k}) ./ rangeVals(k);
        MSE{k} = mean ((meanMatrix{k}-rangeVals(k)).^2 + stdMatrix{k}.^2);
        
        meanNH{k} = mean(WLSnh{k}.range,2);
        stdNH{k} = std(WLSnh{k}.range')';
        MSEnh{k} = mean ((meanNH{k}-rangeVals(k)).^2 + stdNH{k}.^2);

    end
end



% what coefficient gives the lowest bias on average
figure
subplot(1,2,1)
hold on
for k = 1:length(rangeVals)
    ht(k) = plot(WLScoeffVals, BIAS{k}, 'linewidth', 2);
    legendText{k} = ['r = ' num2str(rangeVals(k)) ' km'];
end
xlabel('c [km]')
ylabel('Average bias')
legend(ht, legendText, 'location', 'Northeast')
set(gca, 'yLim', [0 7])
set(gca, 'xLim', [0 12])
text(-0.1,-0.07,'(a)','Units', 'Normalized', 'VerticalAlignment', 'Top')


% what coefficient gives the lowest coefficient of variation on average
subplot(1,2,2)
hold on
for k = 1:length(rangeVals)
    ht(k) = plot(WLScoeffVals, AvgCOV{k}, 'linewidth', 2);
    legendText{k} = ['r = ' num2str(rangeVals(k)) ' km'];
end
xlabel('c [km]')
ylabel('Average coefficient of variation')
set(gca, 'yLim', [0 1])
set(gca, 'xLim', [0 12])
text(-0.1,-0.07,'(b)','Units', 'Normalized', 'VerticalAlignment', 'Top')
Format2x1SubplotFigure
print('-dpdf', [options.figurePath 'WLS_coefficient_study.pdf']); % save the figure to a file




%% plot weighting functional form
figure
h=0:0.5:60;
hold on
plot(h, exp(-h/9), '--b', 'linewidth', 2)
plot(h, exp(-h/5), '-k', 'linewidth', 2)
plot(h, exp(-h/3), '-b', 'linewidth', 2)
legend('c = 9 km', 'c = 5 km', 'c = 3 km')
ylabel('Weight factor')
xlabel('h [km]')
text(-0.1,-0.07,'(b)','Units', 'Normalized', 'VerticalAlignment', 'Top')
FormatFigureBook

%% plots of ranges for each replicate, versus number of stations
yL = [0 120]; % y axis limits
xL = [0 100*ceil(max(numRecs(eventIdx))/100)]; % x axis limits (go to the next round number above max value)

cIdx = [6 10 15 21]; % coeff values to consider
k = 3; % index of range value to consider

figure
for j=1:4 
    subplot(2,2,j)
    hold on
    plot(numRecs(eventIdxTest), WLSsynthetic{k}.range{cIdx(j)}, '.g')
    [xSort,meanSort] = fn_kernel_smooth(numRecs(eventIdxTest), meanMatrix{k}(:,cIdx(j)));
    plot(xSort,meanSort, '-k', 'linewidth', 1)
    [xSort,stdSort] = fn_kernel_smooth(numRecs(eventIdxTest), stdMatrix{k}(:,cIdx(j)));
    plot(xSort,meanSort+stdSort, '-.k', 'linewidth', 1)
    plot(xSort,meanSort-stdSort, '-.k', 'linewidth', 1)

    plot(xL, [1 1]*rangeVals(k), ':k')
    set(gca, 'ylim', [0 yL(2)])
    xlabel(['Number of stations']);
    ylabel('Simulated ranges')
    title(['c = ' num2str(WLScoeffVals(cIdx(j))) ' km'])
    text(-0.1,-0.07,figLabel{j},'Units', 'Normalized', 'VerticalAlignment', 'Top')
end
Format2x2SubplotFigure
print('-dpdf', [options.figurePath 'WLS_coeff_scatter_plots.pdf']); % save the figure to a file

