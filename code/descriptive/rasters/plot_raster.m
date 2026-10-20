function [fh, ts, ax] = plot_raster(cache_name, varargin)
% PLOT_RASTER  Raster of any single-column cache, grouped by moving flies
%   [fh, ts, ax] = plot_raster('freeze')
%   [fh, ts, ax] = plot_raster('speed', 'loom_speed', 25, 'smooth', 30)
%
%   cache_name   name of <cache_name>_cache.mat in the caches folder
%                ('freeze', 'motion', 'speed', 'jumps', 'pixel', ...)
%
%   Top: mean of the cache across the flies of each moving-flies condition.
%   Bottom: one row per focal fly, grouped by condition (coloured bar on the
%   left) and sorted within each group by the time spent in the 'sort_by'
%   cache after frame 18000. Median loom times are drawn on both panels.
%
%   Options (all optional; empty = the preset of the cache, see raster_defaults)
%   loom_speed   loom speed to select the flies [50]
%   genotype     genotype code to select the flies, bouts.genotype (LC11 dataset:
%                2 LC11TNT, 3 EmptyTNT) [all]
%   sort_by     cache used to sort the flies within a group ['freeze']
%   paths        output of path_generator [path_generator]
%   smooth       moving-average window (frames) for the top panel [0]
%   dilate       moving-max window (frames) applied to the raster only, so
%                single-frame events (jumps) survive the downsampling [0]
%   xlim         x limits in frames [16200, last frame]
%   mask_by      cache used to mask the raster, e.g. 'freeze': bins where it is
%                0 are transparent, so the raster is shown only during freezing []
%   cmap, clim   colormap and colour limits of the raster
%   ylabel, ylim, yticks   top panel y label, limits and ticks

opt = inputParser;
addRequired(opt, 'cache_name', @(x) ischar(x) || isstring(x));
addParameter(opt, 'loom_speed', 50);
addParameter(opt, 'genotype', []);
addParameter(opt, 'sort_by', 'freeze');
addParameter(opt, 'paths', []);
addParameter(opt, 'smooth', []);
addParameter(opt, 'dilate', []);
addParameter(opt, 'xlim', []);
addParameter(opt, 'cmap', []);
addParameter(opt, 'clim', []);
addParameter(opt, 'ylabel', []);
addParameter(opt, 'ylim', []);
addParameter(opt, 'yticks', []);
addParameter(opt, 'looms', true);
addParameter(opt, 'mask_by', []);
addParameter(opt, 'block_size', 'string');

parse(opt, cache_name, varargin{:});

p = opt.Results;
cache_name = char(cache_name);

% Fill the options left empty with the preset of the cache
d = raster_defaults(cache_name);
for f = fieldnames(d)'
    if isempty(p.(f{1}))
        p.(f{1}) = d.(f{1});
    end
end

paths = p.paths;
if isempty(paths)
    paths = path_generator;
end

% Load Colors
col_mov = colorcet('I2','N', 5);
col_nloom = cmapper([], 30);

% Load the bouts file to extract the condition of each fly
thresholds = define_thresholds;
bouts = importdata(fullfile(paths.dataset, 'bouts.mat'));
bouts = bouts_formatting(bouts, thresholds);

%  Select loom speed
in_selection = bouts.sloom == p.loom_speed;
if ~isempty(p.genotype)
    in_selection = in_selection & bouts.genotype == p.genotype;
end
selected_flies = unique(bouts.fly(in_selection, :));
n_moving_flies = accumarray(bouts.fly, bouts.moving_flies, [], @unique);
n_moving_flies = n_moving_flies(selected_flies);

data_mat = cache2mat(importdata(fullfile(paths.cache_path, [cache_name '_cache.mat'])), 'selected_flies', selected_flies');
sort_mat = cache2mat(importdata(fullfile(paths.cache_path, [p.sort_by '_cache.mat'])), 'selected_flies', selected_flies');
loom_mat = cache2mat(importdata(fullfile(paths.cache_path, 'loom_cache.mat')), 'selected_flies', selected_flies');

if ~isempty(p.mask_by)
    mask_mat = cache2mat(importdata(fullfile(paths.cache_path, [char(p.mask_by) '_cache.mat'])), 'selected_flies', selected_flies');
end

if size(data_mat, 1) ~= numel(selected_flies)
    error('plot_raster:multiColumn', '%s_cache has more than one column per fly.', cache_name);
end

n_t = size(data_mat, 2);
if isempty(p.xlim)
    p.xlim = [16200, n_t];
end

% Loom Times
loom_ts = diff(loom_mat, [], 2) == 1;
n_flies = size(loom_mat, 1);

loom_times = nan(n_flies, 20);

for f = 1:n_flies
    loom_times(f, :) = find(loom_ts(f, :));
end

median_loom_ts = median(loom_times, 1);

% Construct table
ts = table();
ts.sort_time = sum(sort_mat(:, 18000:end), 2, 'omitnan');
ts.moving_flies = n_moving_flies;
ts.data = data_mat;
if ~isempty(p.mask_by)
    ts.mask = logical(mask_mat);
end
ts.fly = selected_flies;
ts = sortrows(ts, {'moving_flies', 'sort_time'}, 'ascend', 'ComparisonMethod','abs');

groups = unique(ts.moving_flies)';
fps = 60;

% Create figure
fh = figure('color', 'w', 'Position', [100 200 900 800]);
tiledlayout(3, 1, 'TileSpacing', 'compact', 'Padding', 'loose');

% plot
nexttile
hold on

for idx_moving = groups

    avg = mean(ts.data(ts.moving_flies == idx_moving, :), 1, 'omitnan');
    if p.smooth > 1
        avg = movmean(avg, p.smooth);
    end
    plot(1:n_t, avg, 'Color', col_mov(idx_moving + 1,:), 'LineWidth', 1.5)

end

ax(1) = gca;
ax(1).Color = 'none';
apply_generic(ax(1), 'no_x', true, 'xticks', [0, 18000, n_t], 'ytick', p.yticks, 'ylim', p.ylim, 'xlim', p.xlim, 'font_size', 28, 'xpos', 'top', 'line_width', 2.5)
ylabel(p.ylabel, 'FontSize', 28)


% loom lines span the y ticks
yt = ax(1).YTick;
if p.looms
    for i = 1:length(median_loom_ts)
        line([median_loom_ts(i), median_loom_ts(i)], [yt(1), yt(end)], 'Color', col_nloom.vars.nloom(10 + i, :), 'LineWidth', 1.2, 'LineStyle', '-', 'clipping', 'on');
    end
end

nexttile(2,[2,1])
hold on
ax(2) = gca;

raster_data = ts.data;
if p.dilate > 1
    raster_data = movmax(raster_data, p.dilate, 2);
end
im = imagesc(ax(2), raster_data);
if ~isempty(p.mask_by)
    im.AlphaData = ts.mask;     % transparent outside the mask
end
colormap(ax(2), p.cmap);
if ~isempty(p.clim)
    clim(ax(2), p.clim);
end
apply_generic(ax(2), 'xticks', [0, 18000, n_t], 'no_yticks', true, 'ylim', [- 10 n_flies + 10], 'xlim', p.xlim, 'font_size', 28, 'line_width', 2.5)
xticklabels({});

xlabel('Time (min)')
xticklabels(round([0, 18000, n_t] / (60 * fps), 2))
ax(2).XLabel.Position(2) = ax(2).YLim(1) - 0.125 * diff(ax(2).YLim);
ax(2).YLabel.Position(1) = ax(2).YLabel.Position(1) - 1000;

% Add lines for median_loom_ts array
if p.looms
    for i = 1:length(median_loom_ts)
        line([median_loom_ts(i), median_loom_ts(i)], [-30, n_flies + 300], 'Color', col_nloom.vars.nloom(10 + i, :), 'LineWidth', 1.2, 'LineStyle', '-', 'clipping','off');
    end
end

for idx_moving = groups

    if strcmp(p.block_size, 'string')
        fill([ax(2).XLim(1)-200,ax(2).XLim(1)-200,ax(2).XLim(1)-1000,ax(2).XLim(1)-1000], [length(find(ts.moving_flies <= idx_moving)) length(find(ts.moving_flies <= idx_moving - 1)) length(find(ts.moving_flies <= idx_moving - 1))   length(find(ts.moving_flies <= idx_moving))], ...
            col_mov(idx_moving + 1,:), 'EdgeColor','none', 'Clipping', 'off');
        ylabel('Focal Flies')

    elseif strcmp(p.block_size, 'square')
        fill([ax(2).XLim(1)-200,ax(2).XLim(1)-200,ax(2).XLim(1)-3000,ax(2).XLim(1)-3000], [length(find(ts.moving_flies <= idx_moving)) length(find(ts.moving_flies <= idx_moving - 1)) length(find(ts.moving_flies <= idx_moving - 1))   length(find(ts.moving_flies <= idx_moving))], ...
            col_mov(idx_moving + 1,:), 'EdgeColor','none', 'Clipping', 'off');
    else

        ylabel('Focal Flies')

    end
end



linkaxes([ax(1), ax(2)], 'x');

% Fixed layout, identical for every cache, so rasters overlap 100% (e.g. for
% animating freeze vs social motion). The axes leave the tiledlayout, whose
% auto-sizing depends on labels and colorbar, and get fixed positions.
tl = ax(1).Parent;
pos_top = [0.16 0.6757 0.68 0.2478];
pos_ras = [0.16 0.11 0.68 0.5307];
ax(1).Parent = fh;
ax(2).Parent = fh;
delete(tl);
ax(1).Position = pos_top;
ax(2).Position = pos_ras;

% The colorbar floats in the right margin; only shown for 'motion'
ch = colorbar(ax(2), 'Location', 'manual');
ax(2).Position = pos_ras;
ch.Position = [pos_ras(1) + pos_ras(3) + 0.03, pos_ras(2), 0.025, pos_ras(4)];
ch.LineWidth = 1.5;
ch.FontSize = 24;
switch cache_name
    case {'motion'}
        ch.Label.String = 'Social Motion';
        % ch.Label.Position(1) = -1.5;   % move label to the left of the bar
        ax(1).YLim = d.clim;
    otherwise
        ch.Visible = 'off';
end

end

function d = raster_defaults(cache_name)
% Colormap, colour limits and top panel axis of each cache
d.cmap = flipud(cbrewer2('Greys', []));
d.clim = [];
d.ylabel = cache_name;
d.ylim = [];
d.yticks = [];
d.smooth = 0;
d.dilate = 0;

switch cache_name
    case {'freeze', 'freezeframe'}     % 1 = freezing, drawn black
        d.cmap = flipud(gray(256));
        d.clim = [0 1];
        d.ylabel = {'Total Fraction', 'Freezing'};
        d.ylim = [-0.1 1.1];
        d.yticks = [0 1];
    case 'motion'       % social motion
        d.cmap = cbrewer2('Reds', []);
        d.clim = [0 10];
        d.ylabel = {'Social', 'Motion'};
    case 'speed'        % focal speed
        d.cmap = cbrewer2('Greens', []);
        d.clim = [0 25];
        d.ylabel = {'Focal', 'Speed'};
    case 'pixel'        % pixel change
        d.clim = [0 300];
        d.ylabel = {'Pixel', 'Change'};
    case 'jumps'        % rare single-frame events, drawn dark on white
        d.cmap = cbrewer2('Purples', []);
        d.clim = [0 1];
        d.ylabel = {'Fraction', 'Jumping'};
        d.smooth = 120;
        d.dilate = 60;
        
end
end
