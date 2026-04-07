% Double DDM test

clearvars

% Load the table first. We will take advantage of an already existing
% dataset.
threshold_imm = 2; threshold_mob = 2; threshold_pc = 4; id_code = sprintf('imm%d_mob%d_pc%d', threshold_imm, threshold_mob, threshold_pc);
paths = path_generator('folder', fullfile('fitting_freezes', 'le_new'), 'bouts_id', id_code, 'imfirst', false);
bouts = importdata(fullfile(paths.dataset, 'bouts.mat'));
motion_cache = importdata(fullfile(paths.cache_path, 'motion_cache.mat'));

% Here you can change the window of loom-evoked
thresholds = define_thresholds('le_window', struct('le_window_sl', [5 55], 'le_window_fl', [5 55]));
bouts = bouts_formatting(bouts, thresholds);

% Process the data: set thresholds for sm, fs and ln. Set minimum duration.
bouts_proc = data_parser_new(bouts, 'type', 'immobility', 'period', 'loom', 'window', 'le', 'nloom', 2:20);

% Extract the pre-freezing social motion
sm_pre = extract_sm_from_bouts(bouts_proc, 'type', 'onsets', 'output_type', 'mat', 'window', [-60 0]);
sm_freeze = extract_sm_from_bouts(bouts_proc, 'type', 'onlyfreeze', 'output_type', 'mat');
sm_freeze_alternative = extract_sm_from_bouts(bouts_proc, 'type', 'onlyfreeze', 'output_type', 'mat', 'delay', 20);

bouts_proc.smp = mean(sm_pre, 2, 'omitnan');
bouts_proc.sm =  mean(sm_freeze_alternative, 2, 'omitnan');

% Extract Social Motion Array or keep it Stationary
extra.soc_mot_array = extract_sm_from_bouts(bouts_proc, 'type', 'onsets', 'output_type', 'mat', 'window', [0 630]);
% extra.soc_mot_array = bouts_proc.sm  .* ones(size(extra.soc_mot_array));

% Censoring and Truncation point for the fitting
points.censoring = 10.5;
points.truncation = 0.5;

model = 'dddm2';
model_results = run_fitting_newer(bouts_proc, points, model, paths, 'export', true, 'bads_display', 'iter', 'n_bads', 2, 'extra', extra, 'vbmc_exhaustive', false, 'pass_ndt', true);
plot_estimates('results', model_results, 'export', true, 'ylimits', [-2 5])

[fh, ax, ax_inset] = fd_conditions('results', model_results, 'no_y', true, 'vis', true, 'bin_size', 5, 'col', 'gray2', 'censored_inset', true);
overlay_fits(fh, ax, ax_inset, 'results', model_results, 'export', true, 'extra', extra, 'censored_inset', true);
overlay_separate_processes(fh, ax, ax_inset, 'results', model_results, 'export', true, 'extra', extra)


%%

% fh = plot_fit('results', model_results, 'conditions', false, 'export', true, 'bin_size', 1, 'censored_inset', true, 'type', 'continuous');
% fh_conditions = plot_fit('results', model_results, 'conditions', true, 'export', true, 'bin_size', 3 , 'type', 'discrete', 'extra', extra);

plot_estimates('results', model_results, 'export', true, 'ylimits', [-2 10])
[fh, ax, ax_inset] = fd_conditions('results', model_results, 'no_y', true, 'vis', true, 'bin_size', 5, 'col', 'gray2', 'censored_inset', true);
[fh, ax, ax_inset] = fd_conditions_for_cosyne('results', model_results, 'no_y', true, 'vis', true, 'bin_size', 6, 'censored_inset', false,'col', 'gray2', 'export', true);

overlay_fits(fh, ax, ax_inset, 'results', model_results, 'export', true, 'extra', extra, 'censored_inset', true);
overlay_fits_for_cosyne(fh, ax, ax_inset, 'results', model_results, 'export', true, 'extra', extra, 'censored_inset', false, 'col', 'gray2');

overlay_separate_processes(fh, ax, ax_inset, 'results', model_results, 'export', true, 'extra', extra)

%fh = kde_fit('results', model_results, 'conditions', false, 'export', true, 'bin_size', 3, 'censored_inset', true, 'type', 'continuous');