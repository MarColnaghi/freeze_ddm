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

n_moving = 3;

% Load CSV files
list = dir(fullfile('/Users/marcocolnaghi/experimental_data/004--social_ddm/dataset_1/', sprintf('1CS%dNorpA%dLC6ChR_5F-50cm-*.csv', n_moving, 4 - n_moving)));

array_pc = nan(size(list, 1), 36000);

for idx_flies = 1:size(list)
    fly_data = readtable(fullfile(list(idx_flies).folder, list(idx_flies).name));
    imm_frames = fly_data.pixelchange > 3;
    imm_frames = bwareaopen(imm_frames, 1); % Remove small immobile bouts
    imm_frames = ~bwareaopen(~imm_frames, 1); % Remove small moving bouts
    array_pc(idx_flies, :) = imm_frames;
    fprintf('fly number %d \n', idx_flies)
end

%
median_loom_ts = [18023, 19045, 19885, 21026, 22166, 23247, 24147, 24867, 25887, 26728, 27628, 28468, ...
    29309, 29969, 30689, 31770, 32610, 33450, 34470, 35251];



%% Different version

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

%% Create figure
fh = figure('color', 'w', 'Position', [100 200 750 800]);
tl = tiledlayout(3, 1, 'TileSpacing', 'compact', 'Padding', 'loose');

col.n_mov_flies = colorcet('I2','N', 5);
% % plot
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
xlabel('Time (s)')
xticklabels({0 5 10})

for i = 1:length(median_loom_ts)
    line([median_loom_ts(i), median_loom_ts(i)], [0, 1], 'Color', col_nloom.vars.nloom(10 + i, :), 'LineWidth', 1.2, 'LineStyle', '-', 'clipping', 'on');
end

nexttile(2,[2,1])
hold on
ax(2) = gca;

fre_imgsc = imagesc(ax(2), ts.freeze, [0, 1]);
colormap(ax(2), ('gray'));
% apply_generic(ax(2), 'xtick', [0, 18000, size(fr_mat, 2)], 'ytick', [1, size(fr_mat, 1)], 'ylim', [- 10 size(fr_mat, 1) + 10], 'xlim', [16200, size(fr_mat, 2)])
apply_generic(ax(2), 'xticks', [0, 18000, size(fr_mat, 2)], 'no_yticks', true, 'ylim', [- 10 size(fr_mat, 1) + 10], 'xlim', [16200, size(fr_mat, 2)], 'font_size', 32)
xticklabels({});
ylabel('Flies')
% xlabel('Time (s)')

% Add lines for median_loom_ts array
for i = 1:length(median_loom_ts)
    line([median_loom_ts(i), median_loom_ts(i)], [-30, size(fr_mat, 1) + 300], 'Color', col_nloom.vars.nloom(10 + i, :), 'LineWidth', 1.2, 'LineStyle', '-', 'clipping','off');
end

for idx_moving = 0:4
    fill([ax(2).XLim(2)+200,ax(2).XLim(2)+200,ax(2).XLim(2)+1000,ax(2).XLim(2)+1000], [length(find(ts.moving_flies <= idx_moving)) length(find(ts.moving_flies <= idx_moving - 1)) length(find(ts.moving_flies <= idx_moving - 1))   length(find(ts.moving_flies <= idx_moving))], ...
        col.n_mov_flies(idx_moving + 1,:), 'EdgeColor','none', 'Clipping', 'off');
end

linkaxes([ax(1), ax(2)], 'x');

exporter(fh, paths, 'raster_grouped_moving_flies_reoriented.pdf');

% 
% for idx_zooms = [0 18000]
%     xlim([idx_zooms length(ts.freeze)]);
%     if idx_zooms == 0
%         exporter(fh, paths, sprintf('raster_nocb_gen%d.pdf', genotype))
%         exporter(fh, paths, sprintf('raster_nocb_gen%d.png', genotype))
%     else
%         exporter(fh, paths, sprintf('raster_nocb_gen%d_zoom.pdf', genotype))
%         exporter(fh, paths, sprintf('raster_nocb_gen%d_zoom.png', genotype))
%     end
% end

%%

% Select number of moving moving flies

mov_flies = 2;
ts_proc = ts(ts.moving_flies == mov_flies, :);

points.censoring = 10.5;
points.truncation = 0.5;

thresholds = define_thresholds;
thresholds.le_window_fl = [5 40];
thresholds.le_window_sl = [15 50];

bouts = importdata(fullfile(paths.dataset, 'bouts.mat'));
bouts = bouts_formatting(bouts, thresholds);
bouts_proc = data_parser_new(bouts, 'type', 'immobility', 'period', 'loom', 'window', 'le', 'nloom', 1:20, 'min_dur', 30);
bouts_proc = bouts_proc(bouts_proc.sloom == ls & bouts_proc.moving_flies == mov_flies, :);

selected_flies = ts_proc.fly;

% Optional: define time resolution
dt = 0.1;   % adjust to your data (e.g. 1 frame, 0.04 s, etc.)

% Find global max time to size vectors
Tmax = max(bouts_proc.ends);
tvec = 0:dt:Tmax;

% Preallocate output as a struct (one field per fly)
fly_sm = struct;

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
            sm_val = fly_bouts.avg_sm_freeze_norm(b);

            idx = tvec >= onset & tvec <= offset;
            sm_vec(idx) = sm_val;
        end
    end

    % Store output
    fly_sm(f).fly = fly_id;
    fly_sm(f).t = tvec;
    fly_sm(f).avg_sm_vec = sm_vec;
end

n_flies = numel(fly_sm);
n_t     = numel(fly_sm(1).t);

SM_mat = nan(n_flies, n_t);

tvec = fly_sm(1).t;
for f = 1:n_flies
    SM_mat(f, :) = fly_sm(f).avg_sm_vec;
end

fh = figure('color', 'w', 'Position', [100 200 860 290]);
tl = tiledlayout(1, 1, 'TileSpacing', 'compact', 'Padding', 'loose');

col.n_mov_flies = colorcet('I2','N', 5);

nexttile
hold on
ax(2) = gca;

sm_img = imagesc(ax(2), tvec, 1:n_flies, SM_mat);

% Colormap + limits
clim([0 0.8])   % adjust to your avg_sm_freeze_norm range

% % Transparency: hide NaNs
set(sm_img, 'AlphaData', ~isnan(SM_mat));
colormap(ax(2), cbrewer2('Reds', []));

% apply_generic(ax(2), 'xtick', [0, 18000, size(fr_mat, 2)], 'ytick', [1, size(fr_mat, 1)], 'ylim', [- 10 size(fr_mat, 1) + 10], 'xlim', [16200, size(fr_mat, 2)])
apply_generic(ax(2), 'xticks', [0, 18000, size(ts_proc.freeze, 2)], 'no_yticks', true, 'ylim', [-1 size(ts_proc.freeze, 1) + 1], 'xlim', [16200, size(ts_proc.freeze, 2)], 'font_size', 28)
xticklabels([0 5 10]);
ylabel('Flies')
xlabel('Time (s)')
ax(2).XLabel.Position(2) = -10;

%cb1 = colorbar(ax(2)); cb1.Visible = 'off';
cb2 = colorbar(ax(2), 'southoutside', 'FontSize', 18, 'LineWidth', 2, 'TickLength', 0.1, 'TickDirection', 'none');
%text(length(fs_mat), size(fs_mat, 1)/2, 'Avg. Social Motion', 'HorizontalAlignment','center','VerticalAlignment','top','FontSize', 18, 'Rotation', 90)

% Add lines for median_loom_ts array
for i = 1:length(median_loom_ts)
    line([median_loom_ts(i), median_loom_ts(i)], [-4, size(ts_proc.freeze, 1) + 4], 'Color', col_nloom.vars.nloom(10 + i, :), 'LineWidth', 1.5, 'LineStyle', '-', 'clipping', 'off');
end

for idx_moving = 0:4
    fill([ax(2).XLim(2)+200,ax(2).XLim(2)+200,ax(2).XLim(2)+1000,ax(2).XLim(2)+1000], [length(find(ts_proc.moving_flies <= idx_moving)) length(find(ts_proc.moving_flies <= idx_moving - 1)) length(find(ts_proc.moving_flies <= idx_moving - 1))   length(find(ts_proc.moving_flies <= idx_moving))], ...
        col.n_mov_flies(idx_moving + 1,:), 'EdgeColor','none', 'Clipping', 'off');
end

exporter(fh, paths, 'raster_sm.pdf')

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
        bouts_sorted.avg_sm_freeze_norm(b);

end


% QUICK SANITY CHECK
fh = figure('color', 'w', 'Position', [100 200 860 270]);

sm_img = imagesc(freeze_SM_mat);

% Colormap + limits
clim([0 0.8])   % adjust to your avg_sm_freeze_norm range

% % Transparency: hide NaNs
set(sm_img, 'AlphaData', ~isnan(freeze_SM_mat));
colormap(gca, cbrewer2('Reds', []));
set(gca,'color','none')

set(gca, 'YDir', 'normal')
apply_generic(gca, 'no_xticks', true, 'yticks', [1 600])
yticklabels([0 10])

xlabel('Sorted Immobilities')
ylabel('Durations (s)')
%colorbar
cb2 = colorbar(gca, 'horiz', 'Position', [0.6 0.55 0.3 0.05], 'FontSize', 24, 'LineWidth', 2, 'TickLength', 0.3, 'TickDirection', 'none');
cb2.Label.String = 'Average Social Motion';
cb2.Ticks = [0 0.4 0.8];
exporter(fh, paths, 'reordered_bouts_sm.pdf')


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

%%


mov_flies = 2;
ts_proc = ts(ts.moving_flies == mov_flies, :);

points.censoring = 10.5;
points.truncation = 0.5;

thresholds = define_thresholds;
thresholds.le_window_fl = [5 40];
thresholds.le_window_sl = [15 50];

bouts = importdata(fullfile(paths.dataset, 'bouts.mat'));
bouts = bouts_formatting(bouts, thresholds);
bouts_proc = data_parser_new(bouts, 'type', 'immobility', 'period', 'loom', 'window', 'all', 'nloom', 1:20, 'min_dur', 30);
bouts_proc = bouts_proc(bouts_proc.sloom == ls & bouts_proc.moving_flies == mov_flies, :);

selected_flies = ts_proc.fly;

% Optional: define time resolution
dt = 0.1;   % adjust to your data (e.g. 1 frame, 0.04 s, etc.)

% Find global max time to size vectors
Tmax = max(bouts_proc.ends);
tvec = 0:dt:Tmax;

% Preallocate output as a struct (one field per fly)
fly_sm = struct;

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
            sm_val = fly_bouts.avg_sm_freeze_norm(b);

            idx = tvec >= onset & tvec <= offset;
            sm_vec(idx) = sm_val;
        end
    end

    % Store output
    fly_sm(f).fly = fly_id;
    fly_sm(f).t = tvec;
    fly_sm(f).avg_sm_vec = sm_vec;
end

n_flies = numel(fly_sm);
n_t     = numel(fly_sm(1).t);

SM_mat = nan(n_flies, n_t);

tvec = fly_sm(1).t;
for f = 1:n_flies
    SM_mat(f, :) = fly_sm(f).avg_sm_vec;
end

fh = figure('color', 'w', 'Position', [100 200 860 290]);
tl = tiledlayout(1, 1, 'TileSpacing', 'compact', 'Padding', 'loose');

col.n_mov_flies = colorcet('I2','N', 5);

nexttile
hold on
ax(2) = gca;

sm_img = imagesc(ax(2), tvec, 1:n_flies, SM_mat);

% Colormap + limits
clim([0 0.001])   % adjust to your avg_sm_freeze_norm range

% % Transparency: hide NaNs
set(sm_img, 'AlphaData', ~isnan(SM_mat));
colormap(ax(2), cbrewer2('Reds', []));

% apply_generic(ax(2), 'xtick', [0, 18000, size(fr_mat, 2)], 'ytick', [1, size(fr_mat, 1)], 'ylim', [- 10 size(fr_mat, 1) + 10], 'xlim', [16200, size(fr_mat, 2)])
apply_generic(ax(2), 'xticks', [0, 18000, size(ts_proc.freeze, 2)], 'no_yticks', true, 'ylim', [-1 size(ts_proc.freeze, 1) + 1], 'xlim', [16200, size(ts_proc.freeze, 2)], 'font_size', 28)
xticklabels([0 5 10]);
ylabel('Flies')
xlabel('Time (s)')
ax(2).XLabel.Position(2) = -10;

%cb1 = colorbar(ax(2)); cb1.Visible = 'off';
cb2 = colorbar(ax(2), 'southoutside', 'FontSize', 18, 'LineWidth', 2, 'TickLength', 0.1, 'TickDirection', 'none');
%text(length(fs_mat), size(fs_mat, 1)/2, 'Avg. Social Motion', 'HorizontalAlignment','center','VerticalAlignment','top','FontSize', 18, 'Rotation', 90)

% Add lines for median_loom_ts array
for i = 1:length(median_loom_ts)
    line([median_loom_ts(i), median_loom_ts(i)], [-4, size(ts_proc.freeze, 1) + 4], 'Color', col_nloom.vars.nloom(10 + i, :), 'LineWidth', 1.5, 'LineStyle', '-', 'clipping', 'off');
end

for idx_moving = 0:4
    fill([ax(2).XLim(2)+200,ax(2).XLim(2)+200,ax(2).XLim(2)+1000,ax(2).XLim(2)+1000], [length(find(ts_proc.moving_flies <= idx_moving)) length(find(ts_proc.moving_flies <= idx_moving - 1)) length(find(ts_proc.moving_flies <= idx_moving - 1))   length(find(ts_proc.moving_flies <= idx_moving))], ...
        col.n_mov_flies(idx_moving + 1,:), 'EdgeColor','none', 'Clipping', 'off');
end

%exporter(fh, paths, 'raster_sm.pdf')
