
paths = path_generator('folder', 'sims/sanity_checks');

% Select the trial 
idx_trial = 1523;

% Select the parameters of the model
true_bound = 2.1;
beta_drift = 3.2;

% Select the dt for signal upsampling
dt = 1/60;
fps = 1/60;

% Run this now for the timevarying model
points.truncation = 3;
points.censoring = 10.5;
n_T = points.censoring / dt;
time_vector = 0:fps:points.censoring;
time_vector_us = 0:dt:points.censoring;

bouts = importdata(fullfile(paths.dataset, 'bouts.mat'));
bouts_proc = data_parser_new(bouts, 'type', 'immobility', 'period', 'loom', 'window', 'le', 'nloom', 2:20);
bouts_proc.smp = ones(height(bouts_proc), 1);
bouts_proc.intercept = ones(height(bouts_proc), 1);

sm = importdata(fullfile(paths.cache_path, "motion_cache.mat"));
sm_signal = extract_sm_from_bouts(bouts_proc);

sm_signal_us = interp1(time_vector(:), sm_signal', time_vector_us(:), 'nearest');

n_sims = 150000;
mu_tv = sm_signal_us(1:end-1 , idx_trial) * beta_drift;
sim_params_vec = [dt, 1, 0, true_bound, 0];

rt_simulated = nan(n_sims, 1);
rt_mex = nan(n_sims, 1);

for idx_sims = 1:n_sims

    rt_simulated(idx_sims) = drift_diff_new('mu_t', mu_tv, 'theta', true_bound, ...
        'dt', dt, 'T', points.censoring);

    rt_mex(idx_sims) = sim_ddm_seeded(mu_tv, sim_params_vec, 1, idx_sims);

end

rt_mex(isnan(rt_mex)) = 11;
rt_simulated(isnan(rt_simulated)) = 11;
rt_simulated(rt_simulated < points.truncation) = nan;
rt_mex(rt_mex < points.truncation) = nan;

figure('Color', 'w', 'Position', [10 10 610 210])
hold on
histogram(rt_simulated(~isnan(rt_simulated)), 0:1/20:12, ...
    'Normalization', 'pdf', 'FaceColor', 'r', 'EdgeColor', 'none', ...
    'FaceAlpha', 0.3);
histogram(rt_mex(~isnan(rt_mex)), 0:1/20:12, ...
    'Normalization', 'pdf', 'FaceColor', 'b', 'EdgeColor', 'none', ...
    'FaceAlpha', 0.3);

params = param_res(true_bound, 'points', points, 'dt', dt);
[pde_results] = ddm_pdf_from_trace([0, true_bound, 0], mu_tv, params);

plot(pde_results.t(1:end-1), pde_results.pdf, 'k--', 'LineWidth', 1);
stem(pde_results.t(end) + fps, pde_results.survival, 'k')

apply_generic(gca, 'no_y', true)

% Now we will repeat the analysis but with a model with stationary evidence

mu_st =  bouts_proc.sm(idx_trial) * ones(size(mu_tv));

for idx_sims = 1:n_sims

    rt_simulated(idx_sims) = drift_diff_new('mu_t', mu_st, 'theta', true_bound, ...
        'dt', dt, 'T', points.censoring);

    rt_mex(idx_sims) = sim_ddm_seeded(mu_st, sim_params_vec, 1, idx_sims);

end

rt_mex(isnan(rt_mex)) = 11;
rt_simulated(isnan(rt_simulated)) = 11;
rt_simulated(rt_simulated < points.truncation) = nan;
rt_mex(rt_mex < points.truncation) = nan;

figure('Color', 'w', 'Position', [10 10 610 210])
hold on
histogram(rt_simulated(~isnan(rt_simulated)), 0:1/20:12, ...
    'Normalization', 'pdf', 'FaceColor', 'r', 'EdgeColor', 'none', ...
    'FaceAlpha', 0.3);
histogram(rt_mex(~isnan(rt_mex)), 0:1/20:12, ...
    'Normalization', 'pdf', 'FaceColor', 'b', 'EdgeColor', 'none', ...
    'FaceAlpha', 0.3);

params = param_res(true_bound, 'points', points, 'dt', dt);
[pde_results] = ddm_pdf_from_trace([0, true_bound, 0], mu_st, params);

plot(pde_results.t(1:end-1), pde_results.pdf, 'k--', 'LineWidth', 1)
stem(pde_results.t(end) + fps, pde_results.survival, 'k')
pde_results.t = pde_results.t(1:end-1);
trapz(pde_results.t(pde_results.t >= points.truncation & pde_results.t <= points.censoring), pde_results.pdf(pde_results.t >= points.truncation & pde_results.t <= points.censoring)) +  pde_results.survival

apply_generic(gca, 'no_y', true)

[~, f, fd] = nll_fly_ddm_newer([1, true_bound, 0], bouts_proc(idx_trial, :), points, 'model_sddm6', 'iid', 'p', []);
trapz(fd(fd >= points.truncation & fd <= points.censoring), f(fd >= points.truncation & fd <= points.censoring)) +  f(end)

plot(fd, f,  'r.-', 'LineWidth', 1)
