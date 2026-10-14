clearvars

% Load the table first. We will take advantage of an already existing
% dataset.
threshold_imm = 2; threshold_mob = 2; threshold_pc = 4; id_code = sprintf('imm%d_mob%d_pc%d', threshold_imm, threshold_mob, threshold_pc);
paths = path_generator('folder', fullfile('fitting_freezes','le'), 'bouts_id', id_code, 'imfirst', false);

thresholds = define_thresholds;
% thresholds.le_window_fl = [5 40];
% thresholds.le_window_sl = [15 50];

bouts = importdata(fullfile(paths.dataset, 'bouts.mat'));
bouts = bouts_formatting(bouts, thresholds);
bouts_proc = data_parser_new(bouts, 'type', 'immobility', 'period', 'loom', 'window', 'le', 'nloom', 2:20, 'exclude_flies', true);

total_length = 30;

points.censoring = 10.5;
points.truncation = 0.5;

% Mean social motion over the total_length frames before onset, on the sm / 10
% scale. NaN for a bout starting within the first total_length + 1 frames
sm_pre = extract_sm_from_bouts(bouts_proc, 'type', 'onsets', 'output_type', 'mat', 'window', [-total_length -1], 'norm_factor', 10);
bouts_proc.smp = mean(sm_pre, 2);

model_results = run_fitting_newer(bouts_proc, points, 'dddm2', paths, 'export', true, 'bads_display', false, 'n_bads', 5, 'extra', extra, 'vbmc_exhaustive', false, 'pass_ndt', false, 'iid', true);