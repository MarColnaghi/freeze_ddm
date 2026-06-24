%% Script Name:
% rasters
% last update: 21.01

%% Description:
% Generates the rasters

%% Load Data and Preprocess
clear all; close all; clc;

% Colors
col       = cmapper([], 5);
col_nloom = cmapper([], 30);

% Shared data loading (see load_contact_timeseries)
ls   = 50;                       % loom speed
data = load_contact_timeseries(ls);

sm_mat         = data.sm;
fr_mat         = data.fr;
md_mat         = data.md;
selected_flies = data.selected_flies;
n_moving_flies = data.n_moving_flies;

median_loom_ts     = median(data.loom_times, 1);
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

% Contact occupancy: fraction of flies in contact at each frame (any contact,
% regardless of ordering), computed separately for each moving-flies subgroup
% and aligned to the raster timeline below.
in_contact = md_mat < threshold_distance;
conditions = unique(n_moving_flies)';
p_contact  = NaN(numel(conditions), size(md_mat, 2));
for c = 1:numel(conditions)
    p_contact(c, :) = mean(in_contact(n_moving_flies == conditions(c), :), 1);
end
smooth_win = 6;                                   % frames (~0.5 s) rolling average
p_contact  = movmean(p_contact, smooth_win, 2);

% Create figure
fh = figure('color', 'w', 'Position', [100 200 750 850]);
tl = tiledlayout(4, 1, 'TileSpacing', 'compact', 'Padding', 'loose');

col.n_mov_flies = colorcet('I2','N', 5);

xl = [16200, size(md_mat, 2)];

% --- Top tile: population contact fraction over time ---
ax_top = nexttile(1, [1, 1]);
hold on
for c = 1:numel(conditions)
    plot(p_contact(c, :), 'Color', col.n_mov_flies(c, :), 'LineWidth', 1.5)
end
for i = 1:length(median_loom_ts)
    xline(ax_top, median_loom_ts(i), 'Color', col_nloom.vars.nloom(10 + i, :), 'LineWidth', 1, 'HandleVisibility', 'off');
end
apply_generic(ax_top, 'no_xticks', true, 'xlim', xl, 'font_size', 20)
ylabel('Fraction in contact')

% --- Bottom tiles: contact raster ---
ax = nexttile(2, [3, 1]);
hold on

fre_imgsc = imagesc(ax, ts.contact, [0, 1]);
colormap(ax, ('gray'));
apply_generic(ax, 'xticks', [0, 18000, size(fr_mat, 2)], 'no_yticks', true, 'ylim', [- 10 size(fr_mat, 1) + 10], 'xlim', xl, 'font_size', 32)
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

linkaxes([ax_top, ax], 'x');
