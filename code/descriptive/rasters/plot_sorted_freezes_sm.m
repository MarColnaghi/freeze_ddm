function [fh, bouts_sorted] = plot_sorted_freezes_sm(varargin)
% PLOT_SORTED_FREEZES_SM  Loom-evoked freezes sorted by duration, coloured by SM
%   [fh, bouts_sorted] = plot_sorted_freezes_sm('moving_flies', 2)
%
%   One column per freeze (longest on the left), filled up to the observed
%   duration with the average social motion over that observed part. Contact
%   (70 px) and the administrative cap are imposed with impose_contact_threshold.
%   Truncation is applied before and after the contact cut. Freezes cut before the cap (contact or end of the recording) get a black
%   dot at the cut; freezes reaching the cap run to the top of the plot.
%
%   Options (all optional)
%   moving_flies   condition(s), number of moving flies, scalar or vector such as
%                  0:4 [0]. Freezes are grouped by condition, longest first, and a
%                  band below the x axis shows the condition in the plot_raster colours
%   loom_speed     loom speed [50]
%   paths          output of path_generator [path_generator]
%   censoring      administrative cap in s [10.5]
%   truncation     shortest freeze kept in s [0.2]
%   fps            frame rate [60]
%   n_per_group    random subsample of freezes drawn per condition, all of them if
%                  a condition has fewer; [] keeps everything [200]
%   order          'grouped': one block per condition, longest first (as plot_raster);
%                  'mixed': all freezes pooled, longest first;
%                  'interleaved': rank k of every condition side by side, so the
%                  decay of each condition can be compared column by column
%                  ['grouped']
%   seed           seed of the subsampling, so the figure is reproducible [1]

opt = inputParser;
addParameter(opt, 'moving_flies', 0);
addParameter(opt, 'loom_speed', 50);
addParameter(opt, 'paths', []);
addParameter(opt, 'censoring', 10.5);
addParameter(opt, 'truncation', 0.2);
addParameter(opt, 'fps', 60);
addParameter(opt, 'n_per_group', 350);
addParameter(opt, 'seed', 1);
addParameter(opt, 'order', 'grouped', @(x) any(strcmp(x, {'grouped', 'mixed', 'interleaved'})));
parse(opt, varargin{:});
p = opt.Results;

paths = p.paths;
if isempty(paths)
    paths = path_generator;
end

max_len = round(p.censoring * p.fps);
min_len = round(p.truncation * p.fps);

% Loom-evoked freezes of the condition; contact and cap imposed once, so the
% durations are the observed ones and sm the average over the observed part.
% The truncation is then applied once more to the observed time, so a freeze
% touched within truncation s is dropped, as in loom_evoked_bouts
bouts = importdata(fullfile(paths.dataset, 'bouts.mat'));
bouts = bouts_formatting(bouts, define_thresholds);
bouts = data_parser_new(bouts, 'type', 'immobility', 'period', 'loom', 'window', 'le', 'nloom', 1:20, 'min_dur', min_len);
bouts = bouts(bouts.sloom == p.loom_speed & ismember(bouts.moving_flies, p.moving_flies), :);
bouts = impose_contact_threshold(bouts);
bouts = bouts(bouts.durations >= min_len, :);   % truncation again on the observed time

% Same number of freezes per condition: random subsample, reproducible, and the
% global random stream is left as it was
if ~isempty(p.n_per_group)
    stream = RandStream('mt19937ar', 'Seed', p.seed);
    keep = false(height(bouts), 1);
    for g = unique(bouts.moving_flies)'
        idx = find(bouts.moving_flies == g);
        if numel(idx) > p.n_per_group
            idx = idx(randperm(stream, numel(idx), p.n_per_group));
        end
        keep(idx) = true;
    end
    bouts = bouts(keep, :);
end

switch p.order
    case 'grouped'      % block per condition, longest first within it
        bouts_sorted = sortrows(bouts, {'moving_flies', 'durations'}, {'ascend', 'descend'});
    case 'mixed'        % pooled, longest first
        bouts_sorted = sortrows(bouts, 'durations', 'descend');
    case 'interleaved'  % rank within the condition, then condition
        bouts = sortrows(bouts, {'moving_flies', 'durations'}, {'ascend', 'descend'});
        [~, ~, grp] = unique(bouts.moving_flies);
        rank = zeros(height(bouts), 1);
        for g = 1:max(grp)
            rank(grp == g) = 1:nnz(grp == g);
        end
        [~, ord] = sortrows([rank, grp]);
        bouts_sorted = bouts(ord, :);
end
sorted_len = bouts_sorted.durations;
n_freezes = height(bouts_sorted);

freeze_SM_mat = nan(max_len, n_freezes);
len_clamped = min(sorted_len, max_len);
for b = 1:n_freezes
    freeze_SM_mat(1:len_clamped(b), b) = bouts_sorted.sm(b);
end

fh = figure('color', 'w', 'Position', [100 200 1100 400]);
ax = axes('Position', [0.10 0.30 0.85 0.65]);
sm_img = imagesc(freeze_SM_mat);
hold on
clim([0 1])
set(sm_img, 'AlphaData', ~isnan(freeze_SM_mat));
colormap(gca, cbrewer2('Reds', []));

% Freezes cut before the cap: black dot at the cut
cut = bouts_sorted.is_censored & ~bouts_sorted.censored_loom;
plot(find(cut), len_clamped(cut) + 0.5, '|', 'Color', 'k', 'MarkerSize', 2);

set(gca, 'YDir', 'normal')
tick_10s = 10 * p.fps;
apply_generic(gca, 'no_xticks', true, 'yticks', [0.5 tick_10s + 0.5], 'ylims', [0.5 max_len + 0.5], ...
    'font_size', 28, 'line_width', 2.5)
yticklabels([0 10])

xlabel('Sorted Freezes')
ylabel('Durations (s)')
cb = colorbar(gca, 'Location', 'eastoutside', 'FontSize', 24, 'LineWidth', 2, 'TickLength', 0.01, 'TickDirection', 'both');
cb.Label.String = 'Average Social Motion';
cb.Ticks = [0 0.5 1];

% Band below the x axis with the number of moving flies of each freeze, one
% block per run of equal condition, outside the axes (no clipping)
col_mov = colorcet('I2', 'N', 5);
band_y = [-0.02 -0.08] * max_len;
mf = bouts_sorted.moving_flies;
run_start = [1; find(diff(mf) ~= 0) + 1];
run_end = [run_start(2:end) - 1; n_freezes];
n_runs = numel(run_start);
X = [run_start - 0.5, run_end + 0.5, run_end + 0.5, run_start - 0.5]';
Y = repmat([band_y(1); band_y(1); band_y(2); band_y(2)], 1, n_runs);
patch(ax, 'XData', X, 'YData', Y, 'FaceVertexCData', col_mov(mf(run_start) + 1, :), ...
    'FaceColor', 'flat', 'EdgeColor', 'none', 'Clipping', 'off');
xlim(ax, [0.5 n_freezes + 0.5])
ax.XLabel.Position(2) = -0.11 * max_len;

end
