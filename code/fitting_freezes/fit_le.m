% Double DDM test

clearvars

% Load the table first. We will take advantage of an already existing
% dataset.
threshold_imm = 2; threshold_mob = 2; threshold_pc = 4; id_code = sprintf('imm%d_mob%d_pc%d', threshold_imm, threshold_mob, threshold_pc);
paths = path_generator('folder', 'fitting_freezes/le', 'bouts_id', id_code, 'imfirst', false);
bouts = importdata(fullfile(paths.dataset, 'bouts.mat'));
motion_cache = importdata(fullfile(paths.cache_path, 'motion_cache.mat'));

bouts_proc = data_parser_new(bouts, 'type', 'immobility', 'period', 'loom', 'window', 'le');
points.censoring = 10.5;
points.truncation = 0;

kde_estimates = importdata(fullfile('/Users/marcocolnaghi/PhD/freeze_ddm/model_results/fitting_freezes/bsl/kde_spontaneous', id_code, 'kde_estimates_bsl.mat'));
[~, idx] = unique(kde_estimates.Fkde, 'last');
Fkde = kde_estimates.Fkde(idx); xkde = kde_estimates.xkde(idx); fkde = kde_estimates.fkde(idx); 

extra.Fkde = Fkde;
extra.fkde = fkde;
extra.xkde = xkde;

%  Now we added our vector column to the bouts table

chunk_len = points.censoring * 60;
sm_ts = nan(height(bouts_proc), chunk_len);

for idx_trials = 1:height(bouts_proc)

    ons = bouts_proc.onsets(idx_trials);
    off = bouts_proc.ends(idx_trials) - 1;
    sum_motion = motion_cache(bouts_proc.fly(idx_trials));
    sm_raw{idx_trials} = sum_motion(ons:ons + chunk_len - 1) ./ 10;

    sm_ts(idx_trials, :) = sum_motion(ons - 630:ons - 1);

end

soc_mot_array = cell2mat(sm_raw)';
extra.soc_mot_array = soc_mot_array;
results_bsl = importdata('/Users/marcocolnaghi/PhD/freeze_ddm/model_results/fitting_freezes/bsl/exp0/run07/fit_results_exp0.mat');
lambda_est = table2array(results_bsl.estimates_mean);
extra.lambda = lambda_est(~isnan(lambda_est));

%%
model_results = run_fitting_newer(bouts_proc, points, 'expdddm1', paths, 'export', true, 'bads_display', true, 'pass_ndt', false, 'n_bads', 5, 'extra', extra);

%%

fh = plot_fit('results', model_results, 'conditions', false, 'export', true, 'bin_size', 1, 'censored_inset', true, 'type', 'continuous');
fh_conditions = plot_fit('results', model_results, 'conditions', true, 'export', true, 'bin_size', 3 , 'type', 'continuous');
plot_estimates('results', model_results, 'export', true, 'ylimits', [-2 4])

%fh = kde_fit('results', model_results, 'conditions', false, 'export', true, 'bin_size', 3, 'censored_inset', true, 'type', 'continuous');