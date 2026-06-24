clearvars

% Load the table first. We will take advantage of an already existing
% dataset.

%% We need to understand what to do with the censoring.

threshold_imm = 2; threshold_mob = 2; threshold_pc = 4; id_code = sprintf('imm%d_mob%d_pc%d', threshold_imm, threshold_mob, threshold_pc);
paths = path_generator('folder', 'social_motion/proximity', 'bouts_id', id_code, 'imfirst', false);
bouts = importdata(fullfile(paths.dataset, 'bouts.mat'));
motion_cache = importdata(fullfile(paths.cache_path, 'motion_cache.mat'));

% Here you can change the window of loom-evoked
thresholds = define_thresholds('le_window', struct('le_window_sl', [5 55], 'le_window_fl', [5 55]));
% thresholds = define_thresholds;

bouts = bouts_formatting(bouts, thresholds);

% Process the data: set thresholds for sm, fs and ln. Set minimum duration.
trunc_point = 30;

bouts_proc = data_parser_new(bouts, 'type', 'immobility', 'period', 'loom', 'window', 'le', 'nloom', 2:20, 'min_dur', trunc_point);
threshold = 80;
[bouts_proc, contact_mask] = impose_contact_threshold(bouts_proc, 'threshold', threshold);

%%
bouts_proc.ends = bouts_proc.onsets + bouts_proc.censoring_time;
bouts_proc.durations_s = bouts_proc.censoring_time ./ 60;
bouts_proc.durations = bouts_proc.censoring_time;

sm_freeze_full = extract_sm_from_bouts(bouts_proc, 'type', 'onlyfreeze', 'output_type', 'mat','norm_factor', 10, 'cache', 'motion_cache');
window_dur = 15;
sm_freeze_ili = extract_sm_from_bouts(bouts_proc, 'type', 'onsets', 'output_type', 'mat','norm_factor', 10, 'cache', 'motion_cache', 'window', [0 window_dur]);

bouts_proc.sm = mean(sm_freeze_ili, 2, 'omitnan');
% bouts_proc.sm = mean(sm_freeze_full, 2, 'omitnan');

bouts_proc = bouts_proc(bouts_proc.durations > trunc_point, :);
bouts_proc = bouts_proc(bouts_proc.sm < 1.5,:);
type = 'cumulative';
param = {'sm'};

fh = fd_distr_withparam_new('bouts', bouts_proc, 'type', type, 'param', param, 'check_quantiles', true,  'export', true, 'paths', paths);

col = cmapper([], 4);

fh = figure('Position', [100 100 600 450], 'Color', 'w');
t = tiledlayout(3, 1, 'TileSpacing', 'compact', 'Padding', 'loose');
nexttile(1)
quant_list = [0 1/4 1/2 3/4 1];

cens = bouts_proc.contacts ~= 0;          % MATLAB: 1 = CENSORED
edges  = quantile(bouts_proc.sm, quant_list);
sm_bin = discretize(bouts_proc.sm, edges);

% ---- top row: sm distribution + terciles ----
hold on
bin_edges = linspace(min(bouts_proc.sm), max(bouts_proc.sm), 120);
histogram(bouts_proc.sm, 'BinEdges', bin_edges, ...
    'DisplayStyle', 'stairs', 'LineWidth', 1, 'EdgeColor', 'k');
c1 = histcounts(bouts_proc.sm(sm_bin==1), bin_edges);
c2 = histcounts(bouts_proc.sm(sm_bin==2), bin_edges);
c3 = histcounts(bouts_proc.sm(sm_bin==3), bin_edges);
c4 = histcounts(bouts_proc.sm(sm_bin==4), bin_edges);

centers = bin_edges(1:end-1) + diff(bin_edges)/2;

bh = bar(centers, [c1; c2; c3; c4]', 'stacked', 'BarWidth', 1, 'EdgeColor', 'none');
bh(1).FaceColor = col.vars.sm(2,:);
bh(2).FaceColor = col.vars.sm(3,:);
bh(3).FaceColor = col.vars.sm(4,:);
bh(4).FaceColor = col.vars.sm(5,:);

apply_generic(gca, 'no_y', true, 'no_x', true, 'font_size', 20, 'tick_length', 0.015, 'line_width', 2, 'xlim', [0 1.5])

% ---- bottom: censoring-corrected CDFs per tercile ----
nexttile(2, [2, 1])
hold on
for idx_quantiles = 1:4
    sel = sm_bin == idx_quantiles;
    [F, x] = ecdf(bouts_proc.durations_s(sel), 'Censoring', cens(sel));
    plot(x, F, 'LineWidth', 2, 'Color', col.vars.sm(idx_quantiles + 1, :))
end

apply_generic(gca, 'xlim', [0 10.5], 'ylim', [-0.02 1.02], 'font_size', 20, 'tick_length', 0.015, 'line_width', 2)

%% Let's check within the same experimental group

fh = figure('Position', [100 100 1200 450], 'Color', 'w');
t = tiledlayout(3, 5, 'TileSpacing', 'compact', 'Padding', 'loose');

quant_list = [0 1/3 2/3 1];
% Do NOT drop collision bouts — keep them and mark them censored below.
% (removed: bouts_proc = bouts_proc(bouts_proc.contacts == 0, :);)

for idx_mov_flies = 0:4
    bouts_mov = bouts_proc(bouts_proc.moving_flies == idx_mov_flies, :);

    % a collision truncates the bout -> true duration is LONGER -> right-censored
    cens = bouts_mov.contacts ~= 0;          % MATLAB: 1 = CENSORED

    edges  = quantile(bouts_mov.sm, quant_list);
    sm_bin = discretize(bouts_mov.sm, edges);

    % ---- top row: sm distribution + terciles ----
    nexttile(idx_mov_flies + 1)
    hold on
    bin_edges = linspace(min(bouts_proc.sm), max(bouts_proc.sm), 120);
    histogram(bouts_proc.sm, 'BinEdges', bin_edges, ...
        'DisplayStyle', 'stairs', 'LineWidth', 1, 'EdgeColor', 'k');
    c1 = histcounts(bouts_mov.sm(sm_bin==1), bin_edges);
    c2 = histcounts(bouts_mov.sm(sm_bin==2), bin_edges);
    c3 = histcounts(bouts_mov.sm(sm_bin==3), bin_edges);
    
    centers = bin_edges(1:end-1) + diff(bin_edges)/2;

    bh = bar(centers, [c1; c2; c3]', 'stacked', 'BarWidth', 1, 'EdgeColor', 'none');
    bh(1).FaceColor = col.vars.sm(2,:);
    bh(2).FaceColor = col.vars.sm(3,:);
    bh(3).FaceColor = col.vars.sm(4,:);
    apply_generic(gca, 'no_y', true, 'no_x', true, 'font_size', 20, 'tick_length', 0.015, 'line_width', 2, 'xlim', [0 1.5])

    % ---- bottom: censoring-corrected CDFs per tercile ----
    nexttile(idx_mov_flies + 6, [2, 1])
    hold on
    for idx_quantiles = 1:3
        sel = sm_bin == idx_quantiles;
        [F, x] = ecdf(bouts_mov.durations_s(sel), 'Censoring', cens(sel));
        stairs(x, F, 'LineWidth', 2, 'Color', col.vars.sm(idx_quantiles + 1, :))
    end

    apply_generic(gca, 'xlim', [0 10.5], 'ylim', [-0.02 1.02], 'font_size', 20, 'tick_length', 0.015, 'line_width', 2)
end