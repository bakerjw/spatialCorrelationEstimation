function [ ] = fn_magnitude_variability(synthetic, numRecs, magnitude, recIdx, titleText, options, seed)

% Created by Jack Baker, 12/19/2019

figLabel = {'(a)', '(b)', '(c)', '(d)'};

%% plots of ranges for each replicate, versus magnitude
yL = [0 110]; % y axis limits (temporarily hard coded)
xL = [3 8]; % x axis limits (temporarily hard coded)

caseNums = [1 4]; % which fitting cases to plot


figure
for j=1:2
    subplot(2,2,j)
    k = caseNums(j); % set an index for the appropriate fitting case
    hold on
    for i=1:length(recIdx)
        h1 = plot(magnitude(recIdx(i)), synthetic.range{k}(i,seed), '.k');
        yVals(i,1) = synthetic.range{k}(i,seed);
    end
    
    hold on
    b = regress(yVals,[ones(size(yVals)) magnitude(recIdx)]);
    YFIT = b(1) + b(2)*magnitude(recIdx);
    h2 = plot(magnitude(recIdx), YFIT, '-k');

    h3 = plot(xL, [1 1]*options.syntheticRange, ':k');
    set(gca, 'ylim', [0 yL(2)])
    xlabel(['Magnitude']);
    ylabel('Simulated ranges')
    if k==2 
        legend([h1 h2 h3], 'Estimated ranges', 'Regression model', 'Actual range', 'location', 'Northeast')
    end
    title(titleText{k})
    text(-0.1,-0.07,figLabel{j},'Units', 'Normalized', 'VerticalAlignment', 'Top')
end

subplot(2,2,4)
plot(magnitude(recIdx), numRecs, '.k')
xlabel(['Magnitude']);
ylabel(['Number of stations']);
    text(-0.1,-0.07,figLabel{3},'Units', 'Normalized', 'VerticalAlignment', 'Top')

Format2x2SubplotFigure
print('-dpdf', [options.figurePath 'oneRangeReplicateVsMagnitude.pdf']); % save the figure to a file



end

