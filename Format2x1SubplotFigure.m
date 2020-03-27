%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Reformat the currently open figure for consistent display
% last modified by Jack Baker, March 30, 2016
% adapted to work with subplot panels by Peter Stafford, May 13, 2016
% paper size modified by Jack Baker, 12/17/2019
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% target font sizes
labelSize = 9;
axisSize = 9;
legendSize = 8;
annotSize = 7;

% use fixed figure size
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [7 3.25]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperPosition', [0 0 7 3.25]);

% loop over any subplots that exist within
panels = get(gcf,'Children');
for p = 1:length(panels)

    % check on type of child attribute
    if strcmp(get(panels(p),'type'), 'axes')

        % resize axis numbers
        set(panels(p), 'FontSize', axisSize);

        % resize axis labels
        axLabels = get(panels(p),{'XLabel', 'YLabel', 'ZLabel'});
        set([axLabels{:}], 'FontSize', labelSize);

        % resize title text, if there is a title
        titleH = get(panels(p), 'title');
        if ~isempty(titleH)
            if iscell(titleH)
                for i = 1:length(titleH)
                    set(titleH{i}, 'FontSize', labelSize);
                end
            else
                set(titleH, 'FontSize', labelSize);
            end
        end

    elseif strcmp(get(panels(p),'type'), 'legend')
        % resize legend text
        set(panels(p), 'FontSize', legendSize);
    end
end
