paths = path_generator('folder', 'sims/sanity_checks');

% Select the trial 
idx_trial = 3945;

% Select the parameters of the model
true_bound = 2.7;
beta_drift = 4.2;
true_ndt = 0;

% Select the dt for signal upsampling
dt = 1/300;
fps = 1/60;
inizio = 0;
fine = 630;
% Run this now for the timevarying model
points.truncation = 0.5;
points.censoring = 10.5;
n_T = points.censoring / dt;
time_vector = 0:fps:points.censoring;
time_vector_us = 0:dt:points.censoring;

bouts = importdata(fullfile(paths.dataset, 'bouts.mat'));
bouts_proc = data_parser_new(bouts, 'type', 'immobility', 'period', 'loom', 'window', 'le', 'nloom', 2:20);
bouts_proc.smp = ones(height(bouts_proc), 1);
bouts_proc.intercept = ones(height(bouts_proc), 1);

sm = importdata(fullfile(paths.cache_path, "motion_cache.mat"));
sm_signal = extract_sm_from_bouts(bouts_proc, 'type', 'onsets', 'output_type', 'mat', 'window', [inizio fine]);

sm_signal_us = interp1(time_vector(:), sm_signal', time_vector_us(:), 'nearest');

n_sims = 500000;
mu_tv = sm_signal_us(1:end-1 , idx_trial) * beta_drift;
sim_params_vec = [dt, 1, 0, true_bound, 0];

rt_simulated = nan(n_sims, 1);
rt_mex = nan(n_sims, 1);

for idx_sims = 1:n_sims

    rt_mex(idx_sims) = sim_ddm_seeded(mu_tv, sim_params_vec, 1, idx_sims) + true_ndt;

end

rt_mex(isnan(rt_mex)) = 11;
rt_simulated(isnan(rt_simulated)) = 11;
rt_simulated(rt_simulated < points.truncation) = nan;
rt_mex(rt_mex < points.truncation) = nan;