% Double DDM test

clear all

% Load the table first. We will take advantage of an already existing
% dataset.
threshold_imm = 3; threshold_mob = 3; threshold_pc = 4; id_code = sprintf('imm%d_mob%d_pc%d', threshold_imm, threshold_mob, threshold_pc);
paths = path_generator('folder', 'fitting_tests/ed', 'bouts_id', id_code);
load(fullfile(paths.dataset, 'bouts.mat'));
bouts_proc = data_parser_new(bouts, 'type', 'immobility', 'period', 'loom', 'window', 'le');
% bouts_proc = bouts_proc(bouts_proc.nloom < 15, :);

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
    'ground_truth', [0 0], ...
    'link', link_linear ...
    );

% For theta 1
model.theta = struct(...
    'predictors', {{ ...
    struct('name', 'sm'), ...
    struct('name', 'fs'), ...
    struct('name', 'ls'), ...
    struct('name', 'intercept') ...
    }}, ...
    'ground_truth', [-0.25 0.1 0.12 0.25], ...
    'link', link_linear ...
    );

% Non decision time
model.tndt = struct( ...
    'predictors', {{ ...
    struct('name', 'intercept') ...
    }}, ...
    'ground_truth', 0.1, ...
    'link', link_linear ...
    );

[gt, lbl] = get_ground_truth_vector(model);
x = gt(~isnan(gt));
gt_table = array2table(gt, 'VariableNames', lbl);
ncomp_vars = evaluate_model(model, gt_table, y);

% Specify the seed
sim_params.rng = 30;
rng(sim_params.rng);

% General simulation parameters
sim_params.n_trials = height(bouts_proc);
sim_params.dt = 1/60;
sim_params.T = 10.5;
sim_params.time_vector = sim_params.dt:sim_params.dt:sim_params.T;
sim_params.z = 0;

% Simulation settings
sim_params.kde_grid = 0:1/1200:120;
sim_params.eval_trials = sim_params.n_trials;
sim_params.num_sims = 20000;

% Censoring/Truncation
points.censoring = sim_params.T;

% Initialize outputs
rt = table;
rt.ed = nan(sim_params.n_trials, 1);

col = cmapper();

% Load the Motion Cache
motion_cache = importdata(fullfile(paths.cache_path, 'motion_cache.mat'));

values = 7:11;              
pmf = [0.20 0.30 0.35 0.10 0.05];  % Skewed left, sums to 1
fs = 60;
N = sim_params.n_trials;

% Sample
samples = randsample(values, N, true, pmf);
samples_sec = samples / fs;

% Extract Social Motion TimeSeries
chunk_len = length(sim_params.time_vector);
figure
tic
for idx_trials = 1:sim_params.n_trials

    ons = bouts_proc.onsets(idx_trials);

    sum_motion = motion_cache(bouts_proc.fly(idx_trials));
    sm_raw{idx_trials} = sum_motion(ons:ons + chunk_len - 1 ./ 10x);
    sm_chunk = sm_raw{idx_trials};

    mu_tv = gt_table.mu_sm .* sm_chunk;
    mu_st = ncomp_vars.mu(idx_trials);
    theta_s = ncomp_vars.theta(idx_trials);
    tndt_s = ncomp_vars.tndt(idx_trials);

    % Simulate RT from full DDM
    [rt.ed(idx_trials), traj_ed] = extrema_detection_new('mu', mu_st, 'theta', theta_s, ...
        'z', sim_params.z, 'dt', sim_params.dt, 'T', sim_params.T, 'ndt', tndt_s); % samples_sec(idx_trials));

    if idx_trials <= 311 && idx_trials >= 300
        nexttile
        plot_traces(traj_ed, sm_chunk, col, theta_s);
    end
end
toc

rt.ed(rt.ed > points.censoring) = sim_params.T + sim_params.dt;
rt.ed(isnan(rt.ed)) = sim_params.T + sim_params.dt; 
points.censoring = sim_params.T;
points.truncation = 0.3;

fh = figure('color', 'w', 'Position', [100, 100, 600, 300]);
tiledlayout(1, 1, 'TileSpacing', 'compact', 'Padding', 'loose')
nexttile
hold on
histogram(rt.ed, 1/120:1/20:sim_params.T + sim_params.dt * 3, 'EdgeColor', 'none', 'Normalization', 'pdf')
apply_generic(gca)
ylabel('Density'); xlabel('Duration (s)'); ax.YAxis.Visible = 'off';
exporter(fh, paths, 'Durations.pdf')

extra.soc_mot_array = cell2mat(sm_raw)';

%  Now we added our vector column to the bouts table.
bouts_proc.durations_s = rt.ed;
bouts_proc.sm = bouts_proc.avg_sm_freeze_norm;
bouts_proc.smp = bouts_proc.avg_sm_freeze_norm;
bouts_proc.fs = bouts_proc.avg_fs_1s_norm;
bouts_proc.ls = bouts_proc.sloom_norm;
bouts_proc.ln = bouts_proc.nloom_norm;
bouts_proc.intercept = ones(height(y),1);

model_2_fit = 'simed1';
model_results = run_fitting_newer(bouts_proc, points, model_2_fit, paths, 'export', true,  'ground_truth', gt_table, 'bads_display', true, 'pass_ndt', true, 'n_bads', 5, 'extra', extra);
plot_estimates('results', model_results, 'export', true, 'paths', paths)
bouts_proc.sm = bouts_proc.avg_sm_freeze_norm;
bouts_proc.smp = bouts_proc.avg_sm_freeze_norm;
bouts_proc.fs = bouts_proc.avg_fs_freeze_norm;
bouts_proc.ln = bouts_proc.nloom_norm;
bouts_proc.ls = bouts_proc.sloom_norm;
bouts_proc.intercept = ones(height(bouts_proc), 1);
fh = plot_fit('results', model_results, 'conditions', false, 'export', true, 'censored_inset', true, 'bin_size', 1, 'gt', false, 'type', 'discrete');
fh = plot_fit('results', model_results, 'conditions', true, 'export', true, 'bin_size', 5, 'censored_inset', true, 'type', 'discrete');

%%
% 
% n_samples = 300;  % number of samples from q(x)
% p_ndt_accum = zeros(1, 31);
% samples = vbmc_rnd(model_results.vp, n_samples);
% model_func = str2func(strcat('model_', model_2_fit));
% 
% for s = 1:n_samples
%     x_s = samples(s,:);  % depends on how your VI is implemented
%     [~, p_ndt_s, ndt_values] = ndt_loglik(x_s, extra, model_func, points, bouts_proc);
%     p_ndt_accum = p_ndt_accum + p_ndt_s;
% end
% 
% p_ndt_vi = p_ndt_accum / n_samples;
% [~, idx] = max(p_ndt_vi);
% ndt = ndt_values(idx);
% 
% 
% function [log_g, p_ndt_given_data, ndt_values] = ndt_loglik(x, extra, model_func, points, bouts_individual_fly)
% 
% bif = bouts_individual_fly;
% 
% ts = bif.durations_s;
% y = table;
% y.sm = bif.sm;
% y.smp = bif.smp;
% y.fs = bif.fs;
% y.ln = bif.ln;
% y.ls = bif.ls;
% y.intercept = bif.intercept;
% 
% model = model_func();
% [gt, lbl] = get_ground_truth_vector(model);
% lbl = lbl(~isnan(gt));
% gt_table = array2table(x, 'VariableNames', lbl);
% out = evaluate_model(model, gt_table, y);
% 
% if ~isfield(model, 'tndt')
%     out.tndt = zeros(size(out, 1), 1);
% end
% 
% % assess freeze duration categories
% bet = bif.durations_s <= points.censoring;
% abo = bif.durations_s > points.censoring;
% 
% g = zeros(size(ts));
% 
% fs = 60;
% 
% if size(extra.soc_mot_array, 1) == 1
%     out.mu = repmat(extra.soc_mot_array, height(out.theta), 1) .* x(1) .* (1/fs);
% else
%     out.mu = extra.soc_mot_array .* x(1) .* (1/fs);
% end
% 
% n_ndt = 31;
% ndt_values = (0:(n_ndt-1)) / fs; % convert frames to seconds
% ndt_prior = ones(1, n_ndt) / n_ndt;
% 
% [pdf, cdf] = pdf_cdf({'ed'});
% 
% n_trials = length(ts);
% log_liks = zeros(n_trials, n_ndt);
% 
% for i = 1:n_ndt
% 
%     g = zeros(n_trials, 1);      % or eps, depending on your convention
% 
%     current_ndt = ndt_values(i);
% 
%     out.tndt = current_ndt *  ones(length(out.theta), 1);
% 
%     below = bif.durations_s <  out.tndt;
%     bet   = bif.durations_s >= out.tndt & bif.durations_s <= points.censoring;
%     abo   = bif.durations_s >  points.censoring;
% 
%     f = @(ts, inds) pdf.ed(ts, out.theta(inds), out.mu(inds, :), out.tndt(inds), fs);
%     F = @(ts, inds) cdf.ed(ts, out.theta(inds), out.mu(inds, :), out.tndt(inds), fs);
% 
%     if ~isempty(points.truncation)
%         trunc_factor = @(inds) 1 - F(points.truncation, inds);
%     else
%         trunc_factor = @(inds) ones(size(ts(inds)))';
%     end
% 
%     g(bet) = f(ts(bet), bet) ./ trunc_factor(bet);
%     g(abo) = F(points.censoring, abo) ./ trunc_factor(abo);
% 
%     pdf_vals = g;
%     log_liks(:, i) = log(pdf_vals(:) + 1e-10);
% end
% 
% % at the end:
% ll_per_ndt   = sum(log_liks, 1);
% log_prior    = log(ndt_prior);
% log_weighted = log_prior + ll_per_ndt;
% 
% max_ll       = max(log_weighted);
% log_marginal = max_ll + log(sum(exp(log_weighted - max_ll)));
% 
% log_g = log_marginal;
% 
% % posterior over ndt:
% log_p_ndt_given_data = log_weighted - log_marginal;
% p_ndt_given_data = exp(log_p_ndt_given_data);     % [1 x n_ndt]
% % ndt_values you already defined above (0:(n_ndt-1))/fs
% 
% end

% %%
% 
% %%
% estimated_params = table2array(model_results.estimates_mean);
% estimated_params = estimated_params(~isnan(estimated_params));
% fh = figure('color','w','Position',[100,100, 600, 400]);
% histogram(bouts_proc.durations_s, -1/120:1/20:(points.censoring + 2), 'Normalization', 'probability', 'EdgeColor', 'none')
% hold on
% [~, f, fd] = nll_fly_ddm_newer(estimated_params, bouts_proc, points, strcat('model_', model_2_fit), 'iid', 'p', extra);
% plot(mean(reshape(fd(1:end-1), 3, [])), sum(reshape(f(1:end-1), 3, []), 1), 'k--', 'LineWidth', 0.2)
% xlabel('Freeze Duration (s)')
% ylabel('pmf')
% ylim([0 0.005])
% apply_generic(gca)
% 
% 
% %
% % [nll, f, fd] = nll_fly_ddm_newer([60 3.5 0], bouts_proc, points, 'model_ed1', 'iid', 'p', extra);
% 
% %%
% % Step 1: Run VBMC (marginalizing over ndt)
% vp = model_results.vp;
% 
% % Step 2: Sample from the learned posterior q(θ, μ)
% fs = 60;
% n_samples = 100;
% samples = vbmc_rnd(vp, n_samples); % [n_samples x n_params]
% 
% % Step 3: For each sample, compute posterior over ndt
% ndt_values = (0:31) / fs;
% n_ndt = length(ndt_values);
% ndt_posterior = zeros(n_samples, n_ndt);
%         ndt_prior = ones(1, n_ndt) / n_ndt;
% 
% for s = 1:n_samples
%     theta_s = samples(s, 2) * ones(height(bouts_proc), 1);
%     mu_s = samples(s, 1);
%     
%     if size(extra.soc_mot_array, 1) == 1
%         out.mu = repmat(extra.soc_mot_array, height(out.theta), 1) .* mu_s .* (1/fs);
%     else
%         out.mu = extra.soc_mot_array .* mu_s .* (1/fs);
%     end
% 
%     % Compute unnormalized posterior for each ndt value
%     % p(ndt | data, θ, μ) ∝ p(data | θ, μ, ndt) p(ndt)
%     log_unnorm = zeros(1, n_ndt);
%     
%     for i = 1:n_ndt
%         pdf_vals = ed_vectorized_trials(bouts_proc.durations_s - ndt_values(i), ...
%                                               theta_s , out.mu, fs, 'pdf');
%         log_unnorm(i) = sum(log(pdf_vals + 1e-10)) + log(ndt_prior(i));
%     end
%     
%     % Normalize to get p(ndt | data, θ_s, μ_s)
%     log_unnorm = log_unnorm - max(log_unnorm); % stability
%     ndt_posterior(s, :) = exp(log_unnorm) / sum(exp(log_unnorm));
%     
% end
% 
% % Step 4: Average over posterior samples to get p(ndt | data)
% % This is the posterior predictive distribution for ndt
% ndt_posterior_mean = mean(ndt_posterior, 1);
% 
% % Visualize
% figure;
% bar(ndt_values, ndt_posterior_mean);
% xlabel('ndt (seconds)');
% ylabel('Posterior probability');
% title('p(ndt | data)');
% 
% % Point estimate (posterior mode or mean)
% [~, map_idx] = max(ndt_posterior_mean);
% ndt_map = ndt_values(map_idx);
% 
% % Uncertainty quantification
% ndt_samples = zeros(n_samples, 1);
% for s = 1:n_samples
%     ndt_samples(s) = randsample(ndt_values, 1, true, ndt_posterior(s, :));
% end
% ndt_mean = mean(ndt_samples);
% ndt_std = std(ndt_samples);
% ndt_quantiles = quantile(ndt_samples, [0.025, 0.975]);
% 
% fprintf('ndt MAP: %.3f s\n', ndt_map);
% fprintf('ndt Mean: %.3f ± %.3f s\n', ndt_mean, ndt_std);
% fprintf('ndt 95%% CI: [%.3f, %.3f] s\n', ndt_quantiles(1), ndt_quantiles(2));
% 
% 
% 
% 
% %%
% fh = figure('color','w','Position',[100,100, 600, 400]);
% histogram(bouts_proc.durations_s, -1/120:1/20:(points.censoring + 2), 'Normalization', 'probability', 'EdgeColor', 'none')
% hold on
% plot(mean(reshape(fd(1:end-1), 3, [])), sum(reshape(f(1:end-1), 3, []), 1), 'k--', 'LineWidth', 0.2)
% xlabel('Freeze Duration (s)')
% ylabel('pmf')
% ylim([0 0.005])
% apply_generic(gca)
% 
% 
% %%
% [nll, f, fd] = nll_fly_ddm_newer(model_results.starting_position, bouts_proc(1,:), points, 'model_ed1', 'iid', 'p', extra);
% 
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