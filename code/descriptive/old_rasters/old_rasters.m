%% Script Name:
% rasters
% last update: 21.01

%% Description:
% Generates the rasters

%% Load Data and Preprocess
clear all; close all; clc;

% Load Colors
col = cmapper;

% Load Paths
paths = path_generator('folder', 'descriptive/rasters');

% Different version

clearvars

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

% Construct table
ts = table();
ts.freeze_time = sum(fr_mat(:, 18000:end), 2);
ts.moving_flies = n_moving_flies;
ts.freeze = ~fr_mat;
ts.sm = sm_mat;
ts.fly = selected_flies;
ts = sortrows(ts, {'moving_flies', 'freeze_time'}, 'ascend', 'ComparisonMethod','abs');

% Create figure
fh = figure('color', 'w', 'Position', [100 200 750 800]);
tl = tiledlayout(3, 1, 'TileSpacing', 'compact', 'Padding', 'loose');

col.n_mov_flies = colorcet('I2','N', 5);

% plot
nexttile
hold on
 
summed_freezes = zeros(length(unique(ts.moving_flies)), length(ts.freeze));

for idx_moving = unique(n_moving_flies)'
    
    summed_freezes(idx_moving + 1, :) = 1 - sum(ts.freeze(ts.moving_flies == idx_moving,:), 1) ./ length(find(ts.moving_flies == idx_moving));
    plot(1:length(ts.freeze), summed_freezes(idx_moving + 1, :), 'Color', col.n_mov_flies(idx_moving + 1,:), 'LineWidth', 2)

end

ax(1) = gca;
ax(1).Color = 'none';
apply_generic(ax(1), 'no_x', false, 'xticks', [0, 18000, size(fr_mat, 2)], 'ytick', [0 1], 'ylim', [-0.1 1.1], 'xlim', [16200, size(fr_mat, 2)], 'font_size', 32, 'xpos', 'top')
ylabel({'Total Fraction', 'Freezing'}, 'FontSize', 28)
xlabel('Time (min)')
xticklabels({0 5 10})

for i = 1:length(median_loom_ts)
    line([median_loom_ts(i), median_loom_ts(i)], [0, 1], 'Color', col_nloom.vars.nloom(10 + i, :), 'LineWidth', 1.2, 'LineStyle', '-', 'clipping', 'on');
end

nexttile(2,[2,1])
hold on
ax(2) = gca;

fre_imgsc = imagesc(ax(2), ts.freeze, [0, 1]);
colormap(ax(2), ('gray'));
apply_generic(ax(2), 'xticks', [0, 18000, size(fr_mat, 2)], 'no_yticks', true, 'ylim', [- 10 size(fr_mat, 1) + 10], 'xlim', [16200, size(fr_mat, 2)], 'font_size', 32)
xticklabels({});

ylabel('Focal Flies')
set(ax(2) ,'Layer', 'Top')
ax(1).XLabel.Position(2) = 1.25;
ax(2).YLabel.Position(1) = ax(2).YLabel.Position(1) - 1000;

% Add lines for median_loom_ts array
for i = 1:length(median_loom_ts)
    line([median_loom_ts(i), median_loom_ts(i)], [-30, size(fr_mat, 1) + 300], 'Color', col_nloom.vars.nloom(10 + i, :), 'LineWidth', 1.2, 'LineStyle', '-', 'clipping','off');
end

for idx_moving = 0:4
    fill([ax(2).XLim(1)-200,ax(2).XLim(1)-200,ax(2).XLim(1)-1000,ax(2).XLim(1)-1000], [length(find(ts.moving_flies <= idx_moving)) length(find(ts.moving_flies <= idx_moving - 1)) length(find(ts.moving_flies <= idx_moving - 1))   length(find(ts.moving_flies <= idx_moving))], ...
        col.n_mov_flies(idx_moving + 1,:), 'EdgeColor','none', 'Clipping', 'off');
end

linkaxes([ax(1), ax(2)], 'x');

exporter(fh, paths, 'raster_grouped_moving_flies_reoriented.pdf');

%%

% Select number of moving moving flies

mov_flies = 0;
ts_proc = ts(ts.moving_flies == mov_flies, :);

points.censoring = 10.5;
points.truncation = 0.5;

thresholds = define_thresholds;
thresholds.le_window_fl = [5 55];
thresholds.le_window_sl = [5 55];

bouts = importdata(fullfile(paths.dataset, 'bouts.mat'));
bouts = bouts_formatting(bouts, thresholds);
bouts_proc = data_parser_new(bouts, 'type', 'immobility', 'period', 'loom', 'window', 'le', 'nloom', 1:20, 'min_dur', 30);
bouts_proc = bouts_proc(bouts_proc.sloom == ls & bouts_proc.moving_flies == mov_flies, :);

selected_flies = ts_proc.fly;

% Optional: define time resolution
dt = 0.1;   % adjust to your data (e.g. 1 frame, 0.04 s, etc.)

% Find global max time to size vectors
Tmax = max(bouts_proc.ends) - 1;
tvec = 0:dt:Tmax;

% Preallocate output as a struct (one field per fly)
fly_sm = struct;
fs_mat = cache2mat(fs_cache, 'selected_flies', selected_flies');

for f = 1:numel(selected_flies)

    fly_id = selected_flies(f);

    % Get bouts for this fly (may be empty)
    fly_bouts = bouts_proc(bouts_proc.fly == fly_id, :);

    % Initialize vector
    sm_vec = nan(size(tvec));   % or zeros(size(tvec))

    % Only fill if the fly actually has bouts
    if ~isempty(fly_bouts)
        for b = 1:height(fly_bouts)
            onset  = fly_bouts.onsets(b);
            offset = fly_bouts.ends(b);
            sm_val = fly_bouts.sm(b);

            idx = tvec >= onset & tvec <= offset;
            sm_vec(idx) = sm_val;
        end
    end

    % Store output
    fly_sm(f).fly = fly_id;
    fly_sm(f).t = tvec;
    fly_sm(f).avg_sm_vec = sm_vec;
    fly_sm(f).sum_freeze = sum(fly_sm(f).avg_sm_vec > 0);
    fly_sm(f).focal_speed = fs_mat(f,:);
end

% Get the sorted values and the indexing order
[total_freeze_time, sort_idx] = sort([fly_sm.sum_freeze]);

% Reorder the structure fields using that index
% This assumes fly_sm is a 1xN or Nx1 struct array
fly_sm = fly_sm(sort_idx);

n_flies = numel(fly_sm);
n_t     = numel(fly_sm(1).t);

SM_mat = nan(n_flies, n_t);

tvec = fly_sm(1).t;
for f = 1:n_flies
    SM_mat(f, :) = fly_sm(f).avg_sm_vec;
end

% 1. Prepare Focal Speed Matrix
% We need to make sure the focal speed matrix is sorted the same way as fly_sm
FS_mat = nan(n_flies, size(fly_sm(1).focal_speed, 2));
for f = 1:n_flies
    FS_mat(f, :) = fly_sm(f).focal_speed;
end

% Create time vector for focal speed (assuming it spans the full x-axis)
t_fs = linspace(1, size(FS_mat, 2), size(FS_mat, 2));

% 1. Layout Setup
fh = figure('color', 'w', 'Position', [100 200 1100 350]);

% Define a common position for both axes to ensure perfect overlap
axPos = [0.10, 0.25, 0.75, 0.65]; 

% 2. LAYER 1: The Background (Focal Speed - Grayscale)
ax(1) = axes('Position', axPos);
hold(ax(1), 'on');

fs_img = imagesc(ax(1), t_fs, 1:n_flies, FS_mat, [0 25]); 
colormap(ax(1), flipud(cbrewer2('Greys', []))); 
set(fs_img, 'AlphaData', 0.1); % Subtle background

% Formatting the bottom axis
apply_generic(ax(1), 'xticks', [0, 18000, size(FS_mat, 2)], ...
    'no_yticks', false, 'yticks', [1 n_flies], ...
    'ylim', [0.5 n_flies + 0.5], 'xlim', [16200, size(FS_mat, 2)], 'font_size', 28);
xticklabels(ax(1), [0 5 10]);
ylabel(ax(1), 'Focal Flies'); 
xlabel(ax(1), 'Time (min)');

% 3. LAYER 2: The Foreground (Social Motion - Reds)
% Create a new axis on top of the first one
ax(2) = axes('Position', axPos, ...
    'Color', 'none', ...       % Make background transparent
    'XColor', 'none', ...      % Hide redundant X-axis
    'YColor', 'none', ...      % Hide redundant Y-axis
    'HitTest', 'off');         % Let clicks pass through to ax(1) if needed

hold(ax(2), 'on');

sm_img = imagesc(ax(2), tvec, 1:n_flies, SM_mat);
colormap(ax(2), cbrewer2('Reds', 256));
clim(ax(2), [0 1.0]);

% AlphaData: Show Social Motion ONLY where data exists
set(sm_img, 'AlphaData', ~isnan(SM_mat) & (SM_mat > 0));

% Match the limits exactly to ax(1)
ax(2).XLim = ax(1).XLim;
ax(2).YLim = ax(1).YLim;

% 4. Overlays & Metadata
% We draw loom lines and the social indicator on ax(2) so they are on top
for i = 1:length(median_loom_ts)
   line(ax(2), [median_loom_ts(i), median_loom_ts(i)], ax(2).YLim, ...
       'Color', [col_nloom.vars.nloom(10 + i, :), 0.7], 'LineWidth', 1.5, 'Clipping', 'off');
end

% Social condition indicator (the colored bar on the right)
fill(ax(2), [ax(2).XLim(2)+200, ax(2).XLim(2)+1000, ax(2).XLim(2)+1000, ax(2).XLim(2)+200], ...
     [0.5, 0.5, n_flies+0.5, n_flies+0.5], ...
     col.n_mov_flies(mov_flies + 1,:), 'EdgeColor','none', 'Clipping', 'off');

% 5. Link and Export
linkaxes(ax, 'xy'); % Ensure zooming/panning stays in sync
exporter(fh, paths, 'raster_sm_dual_axis.pdf')

%%

% PARAMETERS
fps = 60;
max_len = round(points.censoring * fps);   % 10.5 s -> 630 frames

% USE PRECOMPUTED BOUT DURATIONS
bout_len = bouts_proc.durations;

% SORT BOUTS BY LENGTH
[sorted_len, sort_idx] = sort(bout_len, 'descend');
bouts_sorted = bouts_proc(sort_idx, :);

% PREALLOCATE MATRIX
n_freezes = height(bouts_sorted);
freeze_SM_mat = nan(max_len, n_freezes);

% FILL MATRIX
for b = 1:n_freezes

    % Clamp bout length to censoring
    len_clamped = min(sorted_len(b), max_len);

    % Fill with average social motion for this freeze
    freeze_SM_mat(1:len_clamped, b) = ...
        bouts_sorted.sm(b);

end


% QUICK SANITY CHECK
fh = figure('color', 'w', 'Position', [100 200 1100 290]);

sm_img = imagesc(freeze_SM_mat);

% Colormap + limits
clim([0 1])   % adjust to your avg_sm_freeze_norm range

% % Transparency: hide NaNs
set(sm_img, 'AlphaData', ~isnan(freeze_SM_mat));
colormap(gca, cbrewer2('Reds', []));
set(gca,'color','none')

set(gca, 'YDir', 'normal')
apply_generic(gca, 'no_xticks', true, 'yticks', [1 600], 'font_size', 28)
yticklabels([0 10])

xlabel('Sorted Freezes')
ylabel('Durations (s)')
%colorbar
cb2 = colorbar(gca, 'horiz', 'Position', [0.6 0.55 0.3 0.05], 'FontSize', 24, 'LineWidth', 2, 'TickLength', 0.01, 'TickDirection', 'both');
cb2.Label.String = 'Average Social Motion';
cb2.Ticks = [0 0.5 1];
exporter(fh, paths, 'reordered_bouts_sm_cosyne.pdf')


%%
% exporter(fh, paths, 'raster_grouped_moving_flies_reoriented.pdf');


% Create figure
fh = figure('color', 'w', 'Position', [100 200 700 300]);
ax1 = axes;

% Set colormap
% Sort rows based on non-zero elements
zero_mask = array_pc;
[~, indx] = sort(sum(zero_mask, 2), 'descend');

% Display fs_mat with transparency
sp_imgsc = imagesc(ax1, zero_mask(indx, :), [0, 1]);
set(sp_imgsc, 'AlphaData', zero_mask(indx, :) < 0.5);
%colorbar
hold on;
colormap(ax1, ('gray'));

% Set figure properties
ax1.YColor = 'k'; 
ax1.XColor = 'k';
ax1.FontSize = 32;
yl = ylabel(ax1, 'Flies', 'FontSize', 32);
xl = xlabel(ax1, 'Time (min)', 'FontSize', 32);
ylabel(ax1, '', 'FontSize', 32);
xlabel(ax1, '', 'FontSize', 32);
set(ax1, 'Xtick', [0, median_loom_ts(1), size(zero_mask, 2)], 'XTicklabels', [0, 5, 10], 'YTick', [0, size(zero_mask, 1)]);

set(ax1,'box','off'); 
set(ax1,'linewidth', 3, 'TickDir','out'); % The only other option is 'in'
xlim(ax1, [median_loom_ts(1) - 1000, size(zero_mask, 2)]);
ylim(ax1, [-1 size(zero_mask,1)+2]);

ylabel('Flies'); 
yticks([])
xlabel('Time (s)'); 
apply_generic(ax1, 24)
%     cb1 = colorbar(ax1); cb1.Visible = 'off';
%     cb2 = colorbar(ax, 'eastoutside', 'FontSize', 18, 'LineWidth', 2, 'TickLength', 0.1, 'TickDirection', 'none');
%     text(length(fs_mat), size(fs_mat, 1)/2, 'Avg. Social Motion', 'HorizontalAlignment','center','VerticalAlignment','top','FontSize', 18, 'Rotation', 90)

colormapp = cbrewer2('Purples', 30);
% Add lines for median_loom_ts array

for i = 1:length(median_loom_ts)pow
    line([median_loom_ts(i), median_loom_ts(i)], [-5, size(zero_mask, 1) + 5], 'Color', colormapp(i + 10,:), 'LineWidth', 2.5, 'LineStyle', '-', 'clipping','on');
end

yl.Position(1) = 16300;
xl.Position(2) = 105;
hold off;
uistack(sp_imgsc, 'top');

% Save figure
exporter(fh, paths, sprintf('raster_bw_%d.pdf', n_moving))
exporter(fh, paths, sprintf('raster_bw_%d.eps', n_moving))

