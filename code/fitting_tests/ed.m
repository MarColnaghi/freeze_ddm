% Double DDM test

clear all

% Load the table first. We will take advantage of an already existing
% dataset.
threshold_imm = 3; threshold_mob = 3; threshold_pc = 4; id_code = sprintf('imm%d_mob%d_pc%d', threshold_imm, threshold_mob, threshold_pc);
paths = path_generator('folder', 'fitting_tests/ed', 'bouts_id', id_code);
load(fullfile(paths.dataset, 'bouts.mat'));
bouts_proc = data_parser_new(bouts, 'type', 'immobility', 'period', 'loom', 'window', 'le');
bouts_proc = bouts_proc(bouts_proc.nloom < 15, :);

% Only save useful variables in the table
y = table;
y.sm = bouts_proc.avg_sm_freeze_norm;
y.fs = bouts_proc.avg_fs_1s_norm;
y.ln = bouts_proc.nloom_norm;
y.ls = bouts_proc.sloom_norm;
y.intercept = ones(height(y),1);   % N‑by‑1 column of ones
y.smp = bouts_proc.avg_sm_freeze_norm;
predictors = y.Properties.VariableNames;

ncomp_vars = table();
link_linear = @(x) x;     % log link for bound height
link_logistic = @(x) 1./(1 + exp(-x));     % log link for bound height

% For mu 1
model.mu = struct( ...
    'predictors', {{ ...
    struct('name', 'sm'), ...
    struct('name', 'intercept')}}, ...
    'ground_truth', [60 0], ...
    'link', link_linear ...
    );

% For theta 1
model.theta = struct(...
    'predictors', {{ ...
%     struct('name', 'fs'), ...
%     struct('name', 'ls'), ...
    struct('name', 'intercept') ...
    }}, ...
    'ground_truth', [3.5], ...
    'link', link_linear ...
    );

% % Non decision time
% model.tndt = struct( ...
%     'predictors', {{ ...
%     struct('name', 'intercept') ...
%     }}, ...
%     'ground_truth', 0.155, ...
%     'link', link_linear ...
%     );

[gt, lbl] = get_ground_truth_vector(model);
x = gt(~isnan(gt));
gt_table = array2table(gt, 'VariableNames', lbl);
ncomp_vars = evaluate_model(model, gt_table, y);

% Specify the seed
sim_params.rng = 1;
rng(sim_params.rng);

% General simulation parameters
sim_params.n_trials = height(bouts_proc);
sim_params.dt = 1/60;
sim_params.T = 30;
sim_params.time_vector = sim_params.dt:sim_params.dt:sim_params.T;
sim_params.z = 0;

% Simulation settings
sim_params.kde_grid = 0:1/1200:120;
sim_params.eval_trials = sim_params.n_trials;
sim_params.num_sims = 20000;

% Censoring/Truncation
points.truncation = [];
points.censoring = sim_params.T;

% Initialize outputs
rt = table;
rt.ed = nan(sim_params.n_trials, 1);

col = cmapper();

% Load the Motion Cache
sim_params.motion_cache_path = fullfile(paths.dataset, 'motion_cache.mat');
motion_cache = importdata(sim_params.motion_cache_path);

% Extract Social Motion TimeSeries
chunk_len = length(sim_params.time_vector);
figure
tic
for idx_trials = 1:sim_params.n_trials

    ons = bouts_proc.onsets(idx_trials);

    sum_motion = motion_cache(bouts_proc.fly(idx_trials));
    sm_raw{idx_trials} = sum_motion(ons:ons + chunk_len - 1);
    sm_chunk = sm_raw{idx_trials};

    mu_tv = gt_table.mu_sm .* sm_chunk + gt_table.mu_intercept;
    mu_st = ncomp_vars.mu(idx_trials);
    theta_s = ncomp_vars.theta(idx_trials);
   %  tndt_s = ncomp_vars.tndt(idx_trials);

    % Simulate RT from full DDM
    [rt.ed(idx_trials), traj_ed] = extrema_detection_new('mu_t', mu_tv, 'theta', theta_s, ...
        'z', sim_params.z, 'dt', sim_params.dt, 'T', sim_params.T, 'ndt', 0.15);

    if idx_trials <= 311 && idx_trials >= 300
        nexttile
        plot_traces(traj_ed, sm_chunk, col, theta_s);
    end
end
toc

%%figure;
fh = figure('color', 'w', 'Position', [100, 100, 600, 300]);
tiledlayout(1, 1, 'TileSpacing', 'compact', 'Padding', 'loose')

nexttile
hold on
histogram(rt.ed, 0:1/10:sim_params.T, 'EdgeColor', 'none', 'Normalization', 'pdf')
apply_generic(gca)
xlabel('Duration (s)'); ax.YAxis.Visible = 'off'; 

exporter(fh, paths, 'Durations.pdf')

rt.ed(isnan(rt.ed)) = sim_params.T + 1; 
points.censoring = sim_params.T;
points.truncation = [];

extra.soc_mot_array = cell2mat(sm_raw)';

%  Now we added our vector column to the bouts table.
bouts_proc.durations_s = rt.ed;
bouts_proc.sm = bouts_proc.avg_sm_freeze_norm;
bouts_proc.smp = bouts_proc.avg_sm_freeze_norm;
bouts_proc.fs = bouts_proc.avg_fs_1s_norm;
bouts_proc.ls = bouts_proc.sloom_norm;
bouts_proc.ln = bouts_proc.nloom_norm;
bouts_proc.intercept = ones(height(y),1);

%
model_results = run_fitting_newer(bouts_proc, points, 'ed2', paths, 'export', false, 'extra', extra, 'ground_truth', gt_table);
plot_estimates('results', model_results, 'export', true, 'paths', paths)
bouts_proc.sm = bouts_proc.avg_sm_freeze_norm;
bouts_proc.smp = bouts_proc.avg_sm_freeze_norm;
bouts_proc.fs = bouts_proc.avg_fs_freeze_norm;
bouts_proc.ln = bouts_proc.nloom_norm;
bouts_proc.ls = bouts_proc.sloom_norm;
bouts_proc.intercept = ones(height(bouts_proc), 1);
%
% [nll, f, fd] = nll_fly_ddm_newer([60 3.5 0], bouts_proc, points, 'model_ed1', 'iid', 'p', extra);

%%
% Step 1: Run VBMC (marginalizing over ndt)
vp = model_results.vp;

% Step 2: Sample from the learned posterior q(θ, μ)
fs = 60;
n_samples = 10000;
samples = vbmc_rnd(vp, n_samples); % [n_samples x n_params]

% Step 3: For each sample, compute posterior over ndt
ndt_values = (0:31) / fs;
n_ndt = length(ndt_values);
ndt_posterior = zeros(n_samples, n_ndt);
        ndt_prior = ones(1, n_ndt) / n_ndt;

for s = 1:n_samples
    theta_s = samples(s, 2) * ones(height(bouts_proc), 1);
    mu_s = samples(s, 1);
    
    if size(extra.soc_mot_array, 1) == 1
        out.mu = repmat(extra.soc_mot_array, height(out.theta), 1) .* mu_s .* (1/fs);
    else
        out.mu = extra.soc_mot_array .* mu_s .* (1/fs);
    end

    % Compute unnormalized posterior for each ndt value
    % p(ndt | data, θ, μ) ∝ p(data | θ, μ, ndt) p(ndt)
    log_unnorm = zeros(1, n_ndt);
    
    for i = 1:n_ndt
        pdf_vals = ed_vectorized_trials(bouts_proc.durations_s - ndt_values(i), ...
                                              theta_s , out.mu, fs, 'pdf');
        log_unnorm(i) = sum(log(pdf_vals + 1e-10)) + log(ndt_prior(i));
    end
    
    % Normalize to get p(ndt | data, θ_s, μ_s)
    log_unnorm = log_unnorm - max(log_unnorm); % stability
    ndt_posterior(s, :) = exp(log_unnorm) / sum(exp(log_unnorm));
end

% Step 4: Average over posterior samples to get p(ndt | data)
% This is the posterior predictive distribution for ndt
ndt_posterior_mean = mean(ndt_posterior, 1);

% Visualize
figure;
bar(ndt_values, ndt_posterior_mean);
xlabel('ndt (seconds)');
ylabel('Posterior probability');
title('p(ndt | data)');

% Point estimate (posterior mode or mean)
[~, map_idx] = max(ndt_posterior_mean);
ndt_map = ndt_values(map_idx);

% Uncertainty quantification
ndt_samples = zeros(n_samples, 1);
for s = 1:n_samples
    ndt_samples(s) = randsample(ndt_values, 1, true, ndt_posterior(s, :));
end
ndt_mean = mean(ndt_samples);
ndt_std = std(ndt_samples);
ndt_quantiles = quantile(ndt_samples, [0.025, 0.975]);

fprintf('ndt MAP: %.3f s\n', ndt_map);
fprintf('ndt Mean: %.3f ± %.3f s\n', ndt_mean, ndt_std);
fprintf('ndt 95%% CI: [%.3f, %.3f] s\n', ndt_quantiles(1), ndt_quantiles(2));




%%
fh = figure('color','w','Position',[100,100, 600, 400]);
histogram(bouts_proc.durations_s, -1/120:1/20:(points.censoring + 2), 'Normalization', 'probability', 'EdgeColor', 'none')
hold on
plot(mean(reshape(fd(1:end-1), 3, [])), sum(reshape(f(1:end-1), 3, []), 1), 'k--', 'LineWidth', 0.2)
xlabel('Freeze Duration (s)')
ylabel('pmf')
ylim([0 0.005])
apply_generic(gca)


%%
[nll, f, fd] = nll_fly_ddm_newer(model_results.starting_position, bouts_proc(1,:), points, 'model_ed1', 'iid', 'p', extra);

function fh_traces = plot_traces(traj_ed, sm_chunk, col, theta_s)

ax = gca;
hold on
x_traj = 0:1/60:(length(traj_ed) - 1)/60;
plot(x_traj, traj_ed, 'Color', col.extremadetection, 'LineWidth', 1.5)
plot(x_traj(2:end), sm_chunk, 'k--')
yline(theta_s, 'LineWidth', 2)
hold on
ylim([-theta_s - 1 theta_s + 1])
xlim([0 20])
apply_generic(ax, 20)
xlabel('Time (frames)')
end