clearvars

% Load the table first. We will take advantage of an already existing
% dataset.
threshold_imm = 2; threshold_mob = 2; threshold_pc = 4; id_code = sprintf('imm%d_mob%d_pc%d', threshold_imm, threshold_mob, threshold_pc);
paths = path_generator('folder', fullfile('fitting_freezes','le'), 'bouts_id', id_code, 'imfirst', false);
motion_cache = importdata(fullfile(paths.cache_path, 'motion_cache.mat'));

thresholds = define_thresholds;
% thresholds.le_window_fl = [5 40];
% thresholds.le_window_sl = [15 50];

bouts = importdata(fullfile(paths.dataset, 'bouts.mat'));
bouts = bouts_formatting(bouts, thresholds);
bouts_proc = data_parser_new(bouts, 'type', 'immobility', 'period', 'loom', 'window', 'le', 'nloom', 2:20, 'exclude_flies', true);

total_length = 30;

points.censoring = 10.5;
points.truncation = 0.5;

[mean_sm_before_freeze, mean_sm_during_freeze] = extract_sm_columns(bouts_proc, motion_cache, 'chunk_dur', total_length);
bouts_proc.smp = mean_sm_before_freeze;

model_results = run_fitting_newer(bouts_proc, points, 'dddm2', paths, 'export', true, 'bads_display', false, 'n_bads', 5, 'extra', extra, 'vbmc_exhaustive', false, 'pass_ndt', false, 'iid', true);