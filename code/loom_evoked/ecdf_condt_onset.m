%% Script Name:
% ecdf_condt_onset
% last update: 09.03
%
% ECDF of freeze durations, one curve per onset latency around the loom, one
% column of tiles per loom speed, censored at the first contact or at the end of
% the observation window. The top row is the real looms and the bottom row the
% baseline mock looms, on the same axes and colours, so a latency whose curves
% sit apart from the rest in the top row can be checked against the same latency
% with no loom. The colour bar of the top row greys out the latencies outside the
% loom-evoked window of `thresholds`.
%
%   ecdf_condt_onset                                  % table and window from loom_evoked_data
%   ecdf_condt_onset(bouts, thresholds)
%   ecdf_condt_onset(bouts, thresholds, 'period', 'loom')
%   ecdf_condt_onset(bouts, thresholds, 'bin_size', 5)
%
% bouts is the full table from loom_evoked_data, the bouts plotted are picked
% by loom_evoked_bouts. Pass [] for either input to take it from
% loom_evoked_data.
% period     'both' (default) loom row and baseline row, 'loom' the real looms
%            only, 'bsl' the baseline mock looms only
% bin_size   frames of onset latency pooled into one curve, 1 (default) is one
%            curve per latency. A latency holds about 10 bouts, so the curves
%            of single latencies are noisy. Bins start at the first latency
%            plotted, -10, and a curve is coloured by the middle of its bin

function ecdf_condt_onset(bouts, thresholds, varargin)

opt = inputParser;
addParameter(opt, 'period', 'both', @(x) any(strcmp(x, {'both', 'loom', 'bsl'})));
addParameter(opt, 'bin_size', 1, @(x) isnumeric(x) && isscalar(x) && x >= 1 && x == round(x));
parse(opt, varargin{:});
bin_size = opt.Results.bin_size;

if nargin < 1, bouts = []; end
if nargin < 2, thresholds = []; end

if isempty(bouts) || isempty(thresholds)
    [bouts_default, thresholds_default] = loom_evoked_data;
    if isempty(bouts), bouts = bouts_default; end
    if isempty(thresholds), thresholds = thresholds_default; end
end

if strcmp(opt.Results.period, 'both')
    periods = {'loom', 'bsl'};
else
    periods = {opt.Results.period};
end
period_label = struct('loom', 'loom', 'bsl', 'baseline');

% Number of onsets
n_onsets = 240;
start_position = - 10;

% Latency to row of the colour map, one row per latency. Sequential, so that a
% latency has one colour, unlike the cyclic map this used to have
colour_row = @(latency) latency + abs(start_position) + 1;
cmap = colorcet('L20', 'N', n_onsets + abs(start_position) + 1);
grey = [0.85 0.85 0.85];

fh = figure('color','w','Position',[100, 100, 1200, 400 + 300 * numel(periods)]);
tl = tiledlayout(numel(periods), 2, "TileSpacing","compact","Padding","compact");

for idx_period = 1:numel(periods)
    period = periods{idx_period};

    % The bouts, already contact-censored, and with durations holding the time
    % each was observed for, which is what ecdf takes together with is_censored
    bouts_proc = loom_evoked_bouts(bouts, 'period', period);

    for idx_slooms = [unique(bouts_proc.sloom)]'
        nexttile

        if idx_slooms == 25
            win_field = 'le_window_sl';
        elseif idx_slooms == 50
            win_field = 'le_window_fl';
        end

        bouts_temp = bouts_proc(bouts_proc.sloom == idx_slooms, :);

        hold on
        for bin_first = start_position:bin_size:n_onsets
            bin_last = min(bin_first + bin_size - 1, n_onsets);
            at_onset = bouts_temp.onsets_loomwin >= bin_first & bouts_temp.onsets_loomwin <= bin_last;

            % A bin with no bout, or none that ended on its own, has no curve
            % to draw and ecdf throws on it
            if ~any(at_onset & ~bouts_temp.is_censored)
                continue
            end

            [f, x] = ecdf(bouts_temp.durations(at_onset), 'Censoring', bouts_temp.is_censored(at_onset));
            plot(x, f, 'Color', cmap(colour_row(round((bin_first + bin_last) / 2)), :) , 'LineWidth', 1)
        end

        ax = gca;
        ylim([0 1])
        xlim([0 600])
        apply_generic(ax)
        xlabel('Immobility Duration (frames)')
        ylabel('ecdf')
        tile_label = sprintf('%s, %d cm/s', period_label.(period), idx_slooms);
        if bin_size > 1
            tile_label = sprintf('%s, %d-frame bins', tile_label, bin_size);
        end
        text(590, 0.04, tile_label, 'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom', 'FontSize', 20)

        % The colour bar goes on the top row only, the colours are the same below
        set(ax, 'CLim', [0, size(cmap,1)]);
        bar_cmap = cmap;
        if idx_period == 1
            % Grey out the latencies outside the window. The window is inclusive
            % at both ends, as in bouts_formatting, so its first and last latency
            % keep their colour
            bar_cmap(1:colour_row(thresholds.(win_field)(1)) - 1, :) = repmat(grey, max(colour_row(thresholds.(win_field)(1)) - 1, 0), 1);
            bar_cmap(colour_row(thresholds.(win_field)(2)) + 1:end, :) = repmat(grey, max(size(cmap, 1) - colour_row(thresholds.(win_field)(2)), 0), 1);
        end
        colormap(ax, bar_cmap);

        if idx_period == 1
            cb = colorbar('LineWidth', 2, 'Location', 'northoutside');
            cb.Label.String = 'Immobility Onset';
            cb.FontSize = 20;
            cb.Ticks = [0 30 60 120 240] + abs(start_position);
            cb.TickLabels = {'on', 'off', 60, 120, 240};
        end
    end
end
end
