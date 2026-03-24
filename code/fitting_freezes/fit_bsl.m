% Fit bsl immobilities

clearvars

% Load the table first. We will take advantage of an already existing
% dataset.
threshold_imm = 3; threshold_mob = 3; threshold_pc = 4; id_code = sprintf('imm%d_mob%d_pc%d', threshold_imm, threshold_mob, threshold_pc);
paths = path_generator('folder', 'fitting_freezes/bsl', 'bouts_id', id_code);

motion_cache = importdata(fullfile(paths.cache_path, 'motion_cache.mat'));

load(fullfile(paths.dataset, 'bouts.mat'));
bouts_proc = data_parser_new(bouts, 'type', 'immobility', 'period', 'bsl', 'window', 'le');
points.censoring = 10.5;
points.truncation = min(bouts_proc.durations_s);

total_length = 30;
chunk_len = points.censoring * 60;
sm_pre = nan(height(bouts_proc), total_length);

for idx_trials = 1:height(bouts_proc)

    ons = bouts_proc.onsets(idx_trials);
    off = bouts_proc.ends(idx_trials) - 1;

    sum_motion = motion_cache(bouts_proc.fly(idx_trials));
    sm_pre(idx_trials, :) = sum_motion(ons - total_length:ons - 1) ./ 10;
end

bouts_proc.avg_sm_pre_norm = mean(sm_pre, 2);
model_results = run_fitting_newer(bouts_proc, points, 'exp4', paths, 'export', true, 'extra', []);
plot_estimates('results', model_results, 'export', true, 'ylimits', [-2 10])
[fh, ax, ax_inset] = fd_conditions('results', model_results, 'no_y', true, 'bin_size', 1);
overlay_fits(fh, ax, ax_inset, 'results', model_results, 'export', true)
ylim([0 10])
xlim([0 1])
