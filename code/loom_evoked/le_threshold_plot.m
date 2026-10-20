function [tl, ax] = le_threshold_plot(bouts, thresholds, frame_range, baseline, dur_var, varargin)
% Histogram of immobility onsets around the loom, one tile per loom speed, with
% the loom-evoked window of `thresholds` drawn on top. frame_range selects the
% freezes by duration: a scalar keeps those of at least that many frames, and
% [min max] keeps the ones between the two, both ends included.
%
% baseline is optional: the same selection for the baseline mock looms, drawn as
% an outline over the bars in the same bins and duration range. Both periods
% have the same flies and 20 looms each, so the counts are directly comparable.
%
% dur_var is the bouts column frame_range applies to, 'durations' by default.
% After impose_contact_threshold that column holds the time a bout was observed
% for, so pass 'duration_raw' to sort bouts by how long the freeze really was
%
% By default it opens its own figure. To draw several ranges in one figure, make
% a tiledlayout and pass it as 'parent': the 2 tiles go in a layout nested in
% it, whose position is set by the caller through the returned tl
%
%   t = tiledlayout(1, 3);
%   tl = le_threshold_plot(bouts, thresholds, [30 59], baseline, 'duration_raw', 'parent', t);
%   tl.Layout.Tile = 2;
%
% Name-value options: 'parent' (default [], a new figure), 'font_size' (24),
% 'legend' (true, the loom / baseline legend on the first tile) and 'ylabels'
% (true, the y tick labels). With a parent the x label is left to the caller.

if nargin < 4
    baseline = [];
end

if nargin < 5
    dur_var = 'durations';
end

opt = inputParser;
addParameter(opt, 'parent', [], @(x) isempty(x) || isa(x, 'matlab.graphics.layout.TiledChartLayout'));
addParameter(opt, 'font_size', 24);
addParameter(opt, 'legend', true);
addParameter(opt, 'yaxis', true);
parse(opt, varargin{:});

parent = opt.Results.parent;
font_size = opt.Results.font_size;

if isscalar(frame_range)
    frame_range = [frame_range Inf];
end

if isinf(frame_range(2))
    range_label = sprintf('duration >= %d', frame_range(1));
else
    range_label = sprintf('duration %d to %d', frame_range(1), frame_range(2));
end

e.quantiles = 5;
col = cmapper([], 5);
if isempty(parent)
    fh = figure('color','w','Position',[100,100, 450, 500]);
    tl = tiledlayout(fh, 2, 1, 'TileSpacing', 'tight', 'Padding', 'compact');
else
    tl = tiledlayout(parent, 2, 1, 'TileSpacing', 'tight', 'Padding', 'compact');
end

i = 0; values = [];

for idx_slooms = [unique(bouts.sloom)]'
    i = i + 1;
    nexttile(tl)
    in_range = bouts.(dur_var) >= frame_range(1) & bouts.(dur_var) <= frame_range(2);
    h_handle = histogram(bouts.onsets_loomwin(in_range & bouts.sloom == idx_slooms), -100.5:1:500.5, 'EdgeColor', 'none', 'FaceColor', col.vars.sloom(3*i, :));
    ax(i) = gca;

    if ~isempty(baseline)
        hold on
        in_range_baseline = baseline.(dur_var) >= frame_range(1) & baseline.(dur_var) <= frame_range(2);
        h_baseline = histogram(baseline.onsets_loomwin(in_range_baseline & baseline.sloom == idx_slooms), -100.5:1:500.5, 'DisplayStyle', 'stairs', 'EdgeColor', 'k', 'LineWidth', 1.2);
        if i == 1 && opt.Results.legend
            legend([h_handle h_baseline], {'loom', 'baseline'}, 'Location', 'northeast', 'Box', 'off', 'FontSize', font_size - 4, 'AutoUpdate', 'off');
        end
    end

    if idx_slooms == 25
        xline(thresholds.le_window_sl,'k-.', 'LineWidth', 1.5)
        % ax(i).XAxis.Visible = 'off';
        xticks(sort([0 thresholds.le_window_sl -50 100]))
        xline(mean(thresholds.le_window_sl) - 1.5, 'k--', 'LineWidth', 0.75 )

    elseif idx_slooms == 50
        xline(thresholds.le_window_fl, 'k-.', 'LineWidth', 1.5)
        xticks(sort([0 thresholds.le_window_fl -50 100]))
        xline(mean(thresholds.le_window_fl) - 2.5, 'k--', 'LineWidth', 0.75 )
    end


    
    apply_generic(ax(i), 'font_size', font_size, 'no_y', ~opt.Results.yaxis, 'tick_length', 0.025, 'line_width', 2);
    values = [values; h_handle.Values(:)];
end

linkaxes(ax(:))
%ylim([-2 round(max(values), -1) + 50])
ylim([-5 155])
xlim([-30 120])

if isempty(parent)
    xlabel('Immobility Onsets (frames)', 'FontName', 'Helvetica-Light')
end

