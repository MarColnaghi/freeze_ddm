% Double DDM test

clear all

% Load the table first. We will take advantage of an already existing
% dataset.
threshold_imm = 2; threshold_mob = 2; threshold_pc = 4; id_code = sprintf('imm%d_mob%d_pc%d', threshold_imm, threshold_mob, threshold_pc);
paths = path_generator('folder', 'fitting_freezes/le', 'bouts_id', id_code, 'imfirst', false);
load(fullfile(paths.dataset, 'bouts.mat'));
motion_cache = importdata(fullfile(paths.cache_path, 'motion_cache.mat'));

bouts_proc = data_parser_new(bouts, 'type', 'immobility', 'period', 'loom', 'window', 'le');
points.censoring = 10.5;
points.truncation = 0.35;

%  Now we added our vector column to the bouts table

chunk_len = points.censoring * 60;

for idx_trials = 1:height(bouts_proc)

    ons = bouts_proc.onsets(idx_trials);
    sum_motion = motion_cache(bouts_proc.fly(idx_trials));
    sm_raw{idx_trials} = sum_motion(ons:ons + chunk_len - 1);

end

soc_mot_array = cell2mat(sm_raw)';
extra.soc_mot_array = soc_mot_array;
model_results = run_fitting_newer(bouts_proc, points, 'ded2', paths, 'export', true, 'bads_display', true, 'pass_ndt', true, 'n_bads', 7, 'extra', extra);
model_results.estimates_mean  
plot_fit('results', model_results, 'conditions', false, 'export', true, 'bin_size', 1, 'censored_inset', true)
plot_fit('results', model_results, 'conditions', true, 'export', true, 'bin_size', 10 )
plot_estimates('results', model_results, 'export', true, 'paths', paths)
