% Do posterior analysis of fitted ranges
%
% Created by Yilin Chen, 3/17/2020
% Edited by Jack Baker, last edited 3/25/2020


clear; close all; clc;
load main_data range eventIdx numRecs recIdx recIdxLg magnitude options EQ_name_string EQ_year
load synthetic_data


%% Set which data to consider

k = 2; % use WLS data

rangeData = range(k, :); % use range values estimated from WLS
% rangeData = synthetic.range{k}(:,1)'; % alternative--use range data from simulated data 


nStations = numRecs(eventIdx);



%% distribution estimates

% empirical total variance
sigma_total = std(rangeData);


% standard deviation of fitted range from synthetic data
sigma_1 = std(synthetic.range{k},[],2);

% standard deviation of range prior
sigma_0 = 20; % hard-coded for SA(1s) data from interpretation of prior analysis
%sigma_0 = sqrt(sigma_total^2 - mean(sigma_1.^2));

% mean of range posterior
range_posterior = sigma_0^2 ./ (sigma_0^2 + sigma_1.^2) .* rangeData' + sigma_1.^2 ./ (sigma_0^2 + sigma_1.^2) * mean(rangeData);

% standard deviation of range posterior
sigma_posterior = sqrt(1 ./ (1/sigma_0^2 + 1./sigma_1.^2));


%% what fraction of intervals include the prior mean
fracInInterval = sum(abs(range_posterior-mean(rangeData))./sigma_posterior < 1)/length(range_posterior) 

idxSm = find(numRecs(eventIdx) <= 65);
fracInIntervalSm = sum(abs(range_posterior(idxSm)-mean(rangeData))./sigma_posterior(idxSm) < 1)/length(idxSm) 

idxLg = find(numRecs(eventIdx) > 130);
fracInIntervalSm = sum(abs(range_posterior(idxLg)-mean(rangeData))./sigma_posterior(idxLg) < 1)/length(idxLg) 





%% plot point estimate, and posterior mean and +/- sigma for each earthquake
figure
plot(nStations, rangeData, '.k')
hold on;
errorbar(nStations, range_posterior, sigma_posterior, '+')
plot([40 1000], mean(rangeData)*[1 1], ':k')
xlabel('Number of Stations');
ylabel('Range [km]')
legend('Point estimate of $\hat r$', 'Posterior mean $\pm \sigma$', 'Prior esimate of $\mu_r$', 'location', 'northeast', 'interpreter', 'latex')
set(gca, 'ylim', [-5 120])
set(gca, 'xlim', [0 1000])
set(gca, 'xscale', 'log')
set(gca, 'xtick', [40 100 400 1000])
FormatFigureBook

print('-dpdf', [options.figurePath 'posterior_distributions.pdf']); % save the figure to a file


%% Plot posterior distribution versus point estimate
figure
plot(rangeData, range_posterior, '.b')
hold on
plot([0 125], [0 125], ':k')
plot([0 125], mean(rangeData)*[1 1], '--k')
errorbar(rangeData, range_posterior, sigma_posterior, 'LineStyle','none')
set(gca, 'xlim', [0 125])
xlabel('Point estimate');
ylabel('Posterior distribution')

%% two example posteriors for the two example earthquakes

lineColors{1} = 'b';
lineColors{2} = 'r';

rVals = 0.1:0.1:60;

figure
hold on

for k=1:2

posteriorPDF = normpdf(rVals, range_posterior(options.eventNums(k)), sigma_posterior(options.eventNums(k)) );
h1(k) = plot(rVals, posteriorPDF, 'linewidth', 2, 'color', lineColors{k});
h2(k) = plot(rangeData(options.eventNums(k))*[1 1], [0 0.06], ':', 'linewidth', 1, 'color', lineColors{k});
h3(k) = plot(range_posterior(options.eventNums(k))*[1 1], [0 0.06], '--', 'linewidth', 1, 'color', lineColors{k});

end
legend(h1, 'El Major-Cucapah', 'Yorba Linda')
xlabel('r [km]');
ylabel('Probability density')
FormatFigureBook

print('-dpdf', [options.figurePath 'example_posteriors.pdf']); % save the figure to a file

%% make a summary table of the events (formatted to be a Latex table)

summaryData = [EQ_name_string EQ_year(recIdx) num2cell(eventIdx) num2cell(magnitude(recIdx))  num2cell(numRecs(eventIdx)') num2cell(sigma_1) num2cell(range_posterior)  ]';
sprintf('%s & %.0f & %.0f & %.1f & %.0f & %.1f & %.1f \\\\\n', summaryData{:})

% Header line: \textbf{Event name} & \textbf{Year} & \textbf{Event index} & \textbf{Magnitude}  & \textbf{\# of usable stations} & \textbf{Estimation $\sigma$} & \textbf{Posterior mean $r$ [km]}         \\ 