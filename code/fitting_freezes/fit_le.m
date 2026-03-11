% Double DDM test

clearvars

% Load the table first. We will take advantage of an already existing
% dataset.
threshold_imm = 10; threshold_mob = 10; threshold_pc = 4; id_code = sprintf('imm%d_mob%d_pc%d', threshold_imm, threshold_mob, threshold_pc);
paths = path_generator('folder', fullfile('fitting_freezes', 'le_new'), 'bouts_id', id_code, 'imfirst', false);
bouts = importdata(fullfile(paths.dataset, 'bouts.mat'));
motion_cache = importdata(fullfile(paths.cache_path, 'motion_cache.mat'));

thresholds = define_thresholds;
%thresholds.le_window_fl = [0 60];
%thresholds.le_window_sl = [0 60];

bouts = importdata(fullfile(paths.dataset, 'bouts.mat'));
bouts = bouts_formatting(bouts, thresholds);
bouts_proc = data_parser_new(bouts, 'type', 'immobility', 'period', 'loom', 'window', 'le', 'nloom', 2:20, 'exclude_flies', false);
points.censoring = 10.5;
points.truncation = 0.5;

%  Here we will extract social motion before the freeze onset to feed it to
%  the pmix calculation
sm_pre = extract_sm_from_bouts(bouts_proc, 'type', 'onsets', 'output_type', 'mat', 'window', [-30 -1]);
bouts_proc.smp = mean(sm_pre, 2);

extra = [];
model_results = run_fitting_newer(bouts_proc, points, 'dddm2', paths, 'export', true, 'bads_display', 'none', 'n_bads', 2, 'extra', extra, 'vbmc_exhaustive', false, 'pass_ndt', false);

%%

% fh = plot_fit('results', model_results, 'conditions', false, 'export', true, 'bin_size', 1, 'censored_inset', true, 'type', 'continuous');
% fh_conditions = plot_fit('results', model_results, 'conditions', true, 'export', true, 'bin_size', 3 , 'type', 'discrete', 'extra', extra);

plot_estimates('results', model_results, 'export', true, 'ylimits', [-2 5])
[fh, ax, ax_inset] = fd_conditions('results', model_results, 'no_y', true, 'vis', 'on', 'bin_size', 3);
[fh, ax] = fd_conditions_for_cosyne('results', model_results, 'no_y', true, 'vis', true, 'bin_size', 6, 'censored_inset', false,'col', 'gray', 'export', true);

overlay_fits(fh, ax, ax_inset, 'results', model_results, 'export', true, 'extra', extra, 'censored_inset', false);
overlay_separate_processes(fh, ax, ax_inset, 'results', model_results, 'export', true, 'extra', extra)

%fh = kde_fit('results', model_results, 'conditions', false, 'export', true, 'bin_size', 3, 'censored_inset', true, 'type', 'continuous');