%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Reformat the currently open figure for consistent display
% last modified by Jack Baker, March 30, 2016
% Modified May 19, 2016 by JWB to reduce font sizes
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% target font sizes
labelSize = 9;
axisSize = 9;
legendSize = 8;
annotSize = 7;

% use fixed figure size
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [4.5 3.25]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperPosition', [0 0 4.5 3.25]);

% set font
set(gca, 'FontName', 'Helvetica');

% resize annotation text
textH = findall(findobj(gcf), 'Type', 'text');
set(textH, 'FontSize', annotSize);

% resize axis numbers
set(gca, 'FontSize', axisSize);

% resize axis labels 
axLabels = get(gca,{'XLabel', 'YLabel', 'ZLabel'});
set([axLabels{:}], 'FontSize', labelSize);

% resize legend text
legH = findobj(findobj(gcf), 'tag', 'legend');
set(legH, 'FontSize', legendSize);

% resize title text, if there is a title
titleH = get(gca, 'title');
if ~isempty(titleH)
    if iscell(titleH)
        for i = 1:length(titleH)
            set(titleH{i}, 'FontSize', labelSize);
        end
    else
        set(titleH, 'FontSize', labelSize);
    end
end



