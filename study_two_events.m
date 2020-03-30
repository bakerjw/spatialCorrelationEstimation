% detailed analysis using data from two earthquake events, to illustrate
% some calculations
%
% Created by Jack Baker, 10/8/2019
% Heavily reformatted 12/19/2019
% last edited 3/30/2020

clear; close all; clc;

load singleEventData % load earthquake data
load main_data nPairs h gamma sill range methodName methodNameShort options EQ_name_string % load semivariogram data

options.TStar % output period to screen for checking
EQ_name_string{options.eventNums} % output the earthquake names to the screen





%% plot variograms with one example fit
hPlot = 0:0.5:options.maxR; % finer vector of distances to use for plotting functional predictions

fitIdx = 6; % which fitting method to use for illustration


figure 

k = options.eventNums(1);
subplot(1,2,1) 
plot(h, gamma(k,:), 'ob')
axis([0 options.maxR 0 1.4])
xlabel('h [km]');
ylabel('\gamma(h)');
% title(eventName)
hold on
plot(hPlot, sill(fitIdx,k)*(1-exp(-3*hPlot/range(fitIdx,k))), '-k', 'linewidth', 2)
text(-0.1,-0.07,'(a)','Units', 'Normalized', 'VerticalAlignment', 'Top')
legend('Empirical semivariogram', 'Fitted exponential function', 'location', 'southeast')

k = options.eventNums(2);
subplot(1,2,2) 
plot(h, gamma(k,:), 'ob')
axis([0 options.maxR 0 1.4])
xlabel('h [km]');
ylabel('\gamma(h)');
% title(eventName)
hold on
plot(hPlot, sill(fitIdx,k)*(1-exp(-3*hPlot/range(fitIdx,k))), '-k', 'linewidth', 2)
text(-0.1,-0.07,'(b)','Units', 'Normalized', 'VerticalAlignment', 'Top')
Format2x1SubplotFigure
print('-dpdf', [options.figurePath 'exampleVariogram.pdf']); % save the figure to a file


%% variogram fits from multiple techniques

fitIdxVector = [1, 2, 4, 5]; % which fitting methods to use for illustration

legendText{1} = 'Empirical semivariogram';

k = options.eventNums(1);
figure 
subplot(1,2,1) 
plot(h, gamma(k,:), 'ob')
axis([0 options.maxR 0 1.4])
xlabel('h [km]');
ylabel('\gamma(h)');
hold on
for i = 1:length(fitIdxVector)
    plot(hPlot, sill(fitIdxVector(i),k)*(1-exp(-3*hPlot/range(fitIdxVector(i),k))), '-')
    legendText{i+1} = ['Fitted ' methodNameShort{fitIdxVector(i)} ', $\hat r$ = ' num2str(range(fitIdxVector(i),k),3), ' km'];
end
text(-0.1,-0.07,'(a)','Units', 'Normalized', 'VerticalAlignment', 'Top')
legend(legendText, 'location', 'southeast','interpreter','latex')

k = options.eventNums(2);
subplot(1,2,2) 
plot(h, gamma(k,:), 'ob')
axis([0 options.maxR 0 1.4])
xlabel('h [km]');
ylabel('\gamma(h)');
hold on
for i = 1:length(fitIdxVector)
    plot(hPlot, sill(fitIdxVector(i),k)*(1-exp(-3*hPlot/range(fitIdxVector(i),k))), '-')
    legendText{i+1} = ['Fitted ' methodNameShort{fitIdxVector(i)} ', $\hat r$ = ' num2str(range(fitIdxVector(i),k),3), ' km'];
end
text(-0.1,-0.07,'(b)','Units', 'Normalized', 'VerticalAlignment', 'Top')
legend(legendText, 'location', 'southeast','interpreter','latex')
Format2x1SubplotFigure
print('-dpdf', [options.figurePath 'exampleVariogram4Fits.pdf']); % save the figure to a file

%% plot number of stations at each separation distance 

% compute min and max ratio for the two events
ratio = nPairs(options.eventNums(1),:)./nPairs(options.eventNums(2),:);
minMaxRatio = [min(ratio) max(ratio)] % report ratio of number of stations between the two cases (min and max values)

figure
h1 = plot(h, nPairs, '-c', 'linewidth', 1);
hold on
h2 = plot(h, nPairs(options.eventNums,:), 'linewidth', 2);
xlabel('h [km]');
ylabel('Number of station pairs')
set(gca, 'xlim', [0 options.maxR])
set(gca, 'ylim', [0 700])
legend([h2; h1(1)],  EQ_name_string{options.eventNums}, 'All earthquakes', 'location', 'Northwest')
FormatFigureBook
print('-dpdf', [options.figurePath 'numStationsv2.pdf']); % save the figure to a file



%% data for GMT maps

for k=1:2 % for each event
    % output data files
    data = [longs{options.eventNums(k)}, lats{options.eventNums(k)} resids{options.eventNums(k)}];
    dlmwrite(['data/resids_' num2str(k) '.txt'], data)
    data = [longs{options.eventNums(k)}, lats{options.eventNums(k)} amplitudes{options.eventNums(k)}];
    dlmwrite(['data/SAs_' num2str(k) '.txt'], data)

end




