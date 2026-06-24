%% Script Name:
% rasters
% last update: 21.01

%% Description:
% Generates the rasters

%% Load Data and Preprocess
clear all; close all; clc;

% Load Colors
num_quantiles = 5;
extra.quantiles = num_quantiles;
col = cmapper([], extra.quantiles);
col_nloom = cmapper([], 30);

% Load Paths
paths = path_generator('folder', fullfile('descriptive','rasters'));

% Load timeseries file
sm_cache = importdata(fullfile(paths.cache_path, 'motion_cache.mat'));
pc_cache = importdata(fullfile(paths.cache_path, 'pixel_cache.mat'));
fr_cache = importdata(fullfile(paths.cache_path, 'freeze_cache.mat'));
fs_cache = importdata(fullfile(paths.cache_path, 'speed_cache.mat'));
md_cache = importdata(fullfile(paths.cache_path, 'mindist_cache.mat'));

loom_cache = importdata(fullfile(paths.cache_path, 'loom_cache.mat'));

% Load the bouts file to extract 
threshold_imm = 2; threshold_mob = 2; threshold_pc = 4; id_code = sprintf('imm%d_mob%d_pc%d', threshold_imm, threshold_mob, threshold_pc);
thresholds = define_thresholds;
bouts = importdata(fullfile(paths.dataset, 'bouts.mat'));
bouts = bouts_formatting(bouts, thresholds);

%  Select loom speed
ls = 50;
selected_flies = unique(bouts.fly(bouts.sloom == ls, :));
n_moving_flies = accumarray(bouts.fly, bouts.moving_flies, [], @unique);
n_moving_flies = n_moving_flies(selected_flies);

sm_mat = cache2mat(sm_cache, 'selected_flies', selected_flies');
fs_mat = cache2mat(fs_cache, 'selected_flies', selected_flies');
fr_mat = cache2mat(fr_cache, 'selected_flies', selected_flies');
pc_mat = cache2mat(pc_cache, 'selected_flies', selected_flies');
md_mat = cache2mat(md_cache, 'selected_flies', selected_flies');
loom_mat = cache2mat(loom_cache, 'selected_flies', selected_flies');

% Loom Times
loom_ts = diff(loom_mat, [], 2) == 1;
[r, c] = find(loom_ts);
n_flies = size(loom_mat, 1);

loom_times = nan(n_flies, 20);

for f = 1:n_flies
    loom_times(f, :) = find(loom_ts(f, :));
end

median_loom_ts = median(loom_times, 1);
threshold_distance = 70;

% Construct table
ts = table();
ts.freeze_time = sum(fr_mat(:, 18000:end), 2);
ts.moving_flies = n_moving_flies;
ts.contact = ~(md_mat < threshold_distance);
ts.sum_contact = sum(ts.contact, 2);

ts.sm = sm_mat;
ts.fly = selected_flies;
ts = sortrows(ts, {'moving_flies', 'sum_contact'}, 'ascend', 'ComparisonMethod','abs');

% Create figure
fh = figure('color', 'w', 'Position', [100 200 750 800]);
tl = tiledlayout(2, 1, 'TileSpacing', 'compact', 'Padding', 'loose');

col.n_mov_flies = colorcet('I2','N', 5);

nexttile(1, [2,1])
hold on
ax = gca;

fre_imgsc = imagesc(ax, ts.contact, [0, 1]);
colormap(ax, ('gray'));
apply_generic(ax, 'xticks', [0, 18000, size(fr_mat, 2)], 'no_yticks', true, 'ylim', [- 10 size(fr_mat, 1) + 10], 'xlim', [16200, size(fr_mat, 2)], 'font_size', 32)
xticklabels({});

ylabel('Focal Flies')
set(ax ,'Layer', 'Top')
ax.YLabel.Position(1) = ax.YLabel.Position(1) - 1000;

% Add lines for median_loom_ts array
for i = 1:length(median_loom_ts)
    line([median_loom_ts(i), median_loom_ts(i)], [-30, size(fr_mat, 1) + 300], 'Color', col_nloom.vars.nloom(10 + i, :), 'LineWidth', 1.2, 'LineStyle', '-', 'clipping','off');
end

for idx_moving = 0:4
    fill([ax.XLim(1)-200,ax.XLim(1)-200,ax.XLim(1)-1000,ax.XLim(1)-1000], [length(find(ts.moving_flies <= idx_moving)) length(find(ts.moving_flies <= idx_moving - 1)) length(find(ts.moving_flies <= idx_moving - 1))   length(find(ts.moving_flies <= idx_moving))], ...
        col.n_mov_flies(idx_moving + 1,:), 'EdgeColor','none', 'Clipping', 'off');
end

linkaxes([ax], 'x');
