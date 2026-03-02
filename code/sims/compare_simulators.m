
paths = path_generator('folder', 'sims');
true_bound = 2;
beta_drift = 2.2;
dt = 1/300;
fps = 1/60;

points.truncation = 0;
points.censoring = 10.5;
n_T = points.censoring / dt;
time_vector = 0:fps:points.censoring;
time_vector_us = 0:dt:points.censoring;

bouts = importdata(fullfile(paths.dataset, 'bouts.mat'));
bouts_proc = data_parser_new(bouts, 'type', 'immobility', 'period', 'loom', 'window', 'le', 'nloom', 2:20);
sm = importdata(fullfile(paths.cache_path, "motion_cache.mat"));
sm_signal = extract_sm_from_bouts(bouts_proc);

sm_signal_us = interp1(time_vector(:), sm_signal', time_vector_us(:), 'nearest');

param_res(true_bound, points)

idx_trial = 6786;
n_sims = 50000;
mu_tv = sm_signal_us(1:end-1 , idx_trial) * beta_drift;
sim_params_vec = [dt, 1, 0, true_bound, 0];

for idx_sims = 1:n_sims

    rt_simulated(idx_sims) = drift_diff_new('mu_t', mu_tv, 'theta', true_bound, ...
        'dt', dt, 'T', points.censoring);

    rt_mex(idx_sims) = sim_ddm_seeded(mu_tv, sim_params_vec, 1, idx_sims);


end

figure('Color', 'w')
hold on
histogram(rt_simulated, 0:1/20:points.censoring, 'Normalization', 'pdf')
histogram(rt_mex, 0:1/20:points.censoring, 'Normalization', 'pdf')

params = param_res(true_bound, points);
[pde_results] = ddm_pdf_from_trace([0, true_bound, 0], mu_tv, params);

plot(pde_results.t(1:end-1), pde_results.pdf)