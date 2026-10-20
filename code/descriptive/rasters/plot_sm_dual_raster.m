function [fh, ax, bouts_proc] = plot_sm_dual_raster(ts, varargin)
% PLOT_SM_DUAL_RASTER  Focal speed (grey) + loom-evoked bout social motion (red)
%   [fh, ax, bouts_proc] = plot_sm_dual_raster(ts, 'moving_flies', 2)
%
%   ts           sorted table returned by plot_raster (fly, moving_flies, ...).
%                The rows of the chosen condition keep the order of ts.
%
%   Options (all optional)
%   moving_flies   condition to plot, number of moving flies [0]
%   loom_speed     loom speed used to select the bouts [50]
%   paths          output of path_generator [path_generator]
%   xlim           x limits in frames [16200, last frame]
%   fs_clim        colour limits of the focal speed layer [0 25]
%   sm_clim        colour limits of the social motion layer [0 1]
%
%   Social motion is the average SM of each loom-evoked freeze (data_parser_new:
%   loom period, not frozen at loom start, loom-evoked window, loom 1:20,
%   >= 30 frames), held constant over the bout. Same frame axis and fly order
%   as plot_raster.

opt = inputParser;
addParameter(opt, 'moving_flies', 0);
addParameter(opt, 'loom_speed', 50);
addParameter(opt, 'paths', []);
addParameter(opt, 'xlim', []);
addParameter(opt, 'fs_clim', [0 25]);
addParameter(opt, 'sm_clim', [0 .8]);
parse(opt, varargin{:});
p = opt.Results;

paths = p.paths;
if isempty(paths)
    paths = path_generator;
end

col_mov = colorcet('I2','N', 5);
col_nloom = cmapper([], 30);
fps = 60;

% Caches
fs_cache = importdata(fullfile(paths.cache_path, 'speed_cache.mat'));
loom_mat = cache2mat(importdata(fullfile(paths.cache_path, 'loom_cache.mat')), 'selected_flies', ts.fly');

% Median loom times over all the flies of ts, as in plot_raster
loom_ts = diff(loom_mat, [], 2) == 1;
loom_times = nan(size(loom_mat, 1), 20);
for f = 1:size(loom_mat, 1)
    loom_times(f, :) = find(loom_ts(f, :));
end
median_loom_ts = median(loom_times, 1);

% Flies of the condition, in the order of ts
ts_proc = ts(ts.moving_flies == p.moving_flies, :);
selected_flies = ts_proc.fly;
FS_mat = cache2mat(fs_cache, 'selected_flies', selected_flies');
n_flies = numel(selected_flies);
n_t = size(FS_mat, 2);
tvec = 1:n_t;

if isempty(p.xlim)
    p.xlim = [16200, n_t];
end

% Loom-evoked bouts
thresholds = define_thresholds;
bouts = importdata(fullfile(paths.dataset, 'bouts.mat'));
bouts = bouts_formatting(bouts, thresholds);
bouts_proc = data_parser_new(bouts, 'type', 'immobility', 'period', 'loom', 'window', 'le', 'nloom', 1:20, 'min_dur', 30);
bouts_proc = bouts_proc(bouts_proc.sloom == p.loom_speed & bouts_proc.moving_flies == p.moving_flies, :);

% Bout-averaged SM, held constant over each freeze
SM_mat = nan(n_flies, n_t);
for f = 1:n_flies
    fly_bouts = bouts_proc(bouts_proc.fly == selected_flies(f), :);
    for b = 1:height(fly_bouts)
        idx = tvec >= fly_bouts.onsets(b) & tvec <= fly_bouts.ends(b);
        SM_mat(f, idx) = fly_bouts.sm(b);
    end
end

% Figure: same axis style as plot_raster
fh = figure('color', 'w', 'Position', [100 200 900 340]);
axPos = [0.16, 0.22, 0.68, 0.68];
ylims = [-2 n_flies + 2];
xt = [0, 18000, n_t];

% Layer 1: focal speed
ax(1) = axes('Position', axPos);
hold(ax(1), 'on');
fs_img = imagesc(ax(1), tvec, 1:n_flies, FS_mat, p.fs_clim);
colormap(ax(1), flipud(cbrewer2('Greys', [])));
set(fs_img, 'AlphaData', 0.1);
apply_generic(ax(1), 'xticks', xt, 'no_yticks', true, 'ylims', ylims, ...
    'xlims', p.xlim, 'font_size', 28, 'line_width', 2.5);
xticklabels(ax(1), round(xt / (60 * fps), 2));
xlabel(ax(1), 'Time (min)');
ylabel(ax(1), 'Focal Flies');
ax(1).YLabel.Position(1) = ax(1).YLabel.Position(1) - 1000;

% Layer 2: social motion, transparent axes on top
ax(2) = axes('Position', axPos, 'Color', 'none', 'XColor', 'none', ...
    'YColor', 'none', 'HitTest', 'off');
hold(ax(2), 'on');
sm_img = imagesc(ax(2), tvec, 1:n_flies, SM_mat);
colormap(ax(2), cbrewer2('Reds', 256));
clim(ax(2), p.sm_clim);
set(sm_img, 'AlphaData', ~isnan(SM_mat) & (SM_mat > 0));
ax(2).XLim = ax(1).XLim;
ax(2).YLim = ax(1).YLim;

% Median loom times
for i = 1:numel(median_loom_ts)
    line(ax(2), [median_loom_ts(i), median_loom_ts(i)], [-5, n_flies + 5], ...
        'Color', col_nloom.vars.nloom(10 + i, :), 'LineWidth', 1.2, 'Clipping', 'off');
end

% Condition bar on the left, same size as plot_raster
fill(ax(2), [ax(2).XLim(1)-200, ax(2).XLim(1)-200, ax(2).XLim(1)-1000, ax(2).XLim(1)-1000], ...
    [0.5, n_flies + 0.5, n_flies + 0.5, 0.5], col_mov(p.moving_flies + 1, :), ...
    'EdgeColor', 'none', 'Clipping', 'off');

linkaxes(ax, 'xy');

% Colorbar of the SM layer, floating in the right margin
ch = colorbar(ax(2), 'Location', 'manual');
ax(1).Position = axPos;     % colorbar() can resize its parent axes
ax(2).Position = axPos;
ch.Position = [axPos(1) + axPos(3) + 0.03, axPos(2), 0.025, axPos(4)];
ch.LineWidth = 1.5;
ch.FontSize = 24;
ch.Ticks = p.sm_clim;
ch.Label.String = 'Avg. Social Motion';
end
