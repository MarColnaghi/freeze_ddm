% Double DDM test

clear all

% Load the table first. We will take advantage of an already existing
% dataset.
threshold_imm = 2; threshold_mob = 2; threshold_pc = 4; id_code = sprintf('imm%d_mob%d_pc%d', threshold_imm, threshold_mob, threshold_pc);
paths = path_generator('folder', fullfile('fitting_tests','dddm'), 'bouts_id', id_code);
load(fullfile(paths.dataset, 'bouts.mat'));
bouts_proc = data_parser_new(bouts, 'type', 'immobility', 'period', 'loom', 'window', 'le');
motion_cache = importdata(fullfile(paths.cache_path, 'motion_cache.mat'));

% % Only save useful variables in the table
% y = table;
% y.sm = bouts_proc.avg_sm_freeze_norm;
% y.fs = bouts_proc.avg_fs_1s_norm;
% y.ln = bouts_proc.nloom_norm;
% y.ls = bouts_proc.sloom_norm;
% y.intercept = ones(height(y),1);   % N‑by‑1 column of ones
% predictors = y.Properties.VariableNames;

ncomp_vars = table();
link_linear = @(x) x;     % log link for bound height
link_logistic = @(x) 1./(1 + exp(-x));     % log link for bound height

% For mu 1
model.mu1 = struct( ...
    'predictors', {{ ...
    struct('name', 'sm') ...
    }}, ...
    'ground_truth', [1.3], ...
    'link', link_linear ...
    );

% For theta 1
model.theta1 = struct( ...
    'predictors', {{struct('name', 'intercept')
    }}, ...
    'ground_truth', 0.6, ...
    'link', link_linear ...
    );

% For mu 2
model.mu2 = struct( ...
    'predictors', {{ ...
    struct('name', 'sm') ...
    }}, ...
    'ground_truth', [1.4], ...
    'link', link_linear ...
    );

% For theta
model.theta2 = struct( ...
    'predictors', {{ ...
    struct('name', 'ls'), ...
    struct('name', 'intercept') ...
    }}, ...
    'ground_truth', [0.7 2.2], ...
    'link', link_linear ...
    );

model.pmix = struct( ...
    'predictors', {{ ...
    struct('name', 'fs') ...
    struct('name', 'ls') ...
    }}, ...
    'ground_truth', [-1.0 1.0], ...
    'link', link_logistic ...
    );

% Non decision time
model.tndt = struct( ...
    'predictors', {{ ...
    struct('name', 'intercept') ...
    }}, ...
    'ground_truth', 0.3, ...
    'link', link_linear ...
    );

% Specify the seed
sim_params.rng = 213;
rng(sim_params.rng);

% General simulation parameters
sim_params.n_trials = height(bouts_proc);
sim_params.dt = 1/60;
sim_params.T = 10.5;
sim_params.time_vector = 0:sim_params.dt:sim_params.T;
sim_params.z = 0;

points.censoring = sim_params.T;
points.truncation = 0;

% Simulation settings
sim_params.kde_grid = 0:1/1200:120;
sim_params.eval_trials = sim_params.n_trials;
sim_params.num_sims = 20000;

% Initialize outputs
rt = table;
rt.st = nan(sim_params.n_trials, 1);
rt.ig = nan(sim_params.n_trials, 1);
trial_type = nan(sim_params.n_trials, 1);

% Extract Social Motion TimeSeries
total_length = 30;
chunk_len = points.censoring * 60;
sm_pre = nan(height(bouts_proc), total_length);

for idx_trials = 1:height(bouts_proc)

    ons = bouts_proc.onsets(idx_trials);
    off = bouts_proc.ends(idx_trials) - 1;
    sum_motion = motion_cache(bouts_proc.fly(idx_trials));
    sm_during{idx_trials} = sum_motion(ons:ons + chunk_len) ./ 10;

    sum_motion = motion_cache(bouts_proc.fly(idx_trials));

    sm_pre(idx_trials, :) = sum_motion(ons - total_length:ons - 1) ./ 10;
end

% Add the social motion
soc_mot_array = cell2mat(sm_during)';
extra.soc_mot_array = soc_mot_array;
y.smp = mean(sm_pre, 2);

bouts_proc.intercept = ones(height(bouts_proc), 1);
[gt, lbl] = get_ground_truth_vector(model);
x = gt(~isnan(gt));
gt_table = array2table(gt, 'VariableNames', lbl);
ncomp_vars = evaluate_model(model, gt_table, bouts_proc);

tic
for idx_trials = 1:sim_params.n_trials

    % Determine model 1 or 2 based on pmix
    if rand < ncomp_vars.pmix(idx_trials)
        trial_type(idx_trials) = 1;
        mu_s = ncomp_vars.mu1(idx_trials);
        theta_s = ncomp_vars.theta1(idx_trials);

    else
        trial_type(idx_trials) = 2;
        mu_s = ncomp_vars.mu2(idx_trials);
        theta_s = ncomp_vars.theta2(idx_trials);

    end
    tndt_s = ncomp_vars.tndt(idx_trials);

    % Simulate RT from full DDM
    [rt.st(idx_trials), traj_st] = drift_diff_new('mu', mu_s, 'theta', theta_s, ...
        'z', sim_params.z, 'dt', sim_params.dt, 'T', sim_params.T, 'ndt', tndt_s);

%     % ---- Apply drift scale ----
%     sim_signal_coarse = ones(length(sim_params.dt)) * mu_s;
% 
%     % ---- Interpolate to FINE grid ----
%     fine_signal = interp1(t_coarse(:), sim_signal_coarse(:), ...
%                           t_fine(:), 'nearest');
%     fine_signal(isnan(fine_signal)) = sim_signal_coarse(end);
% 
%     % ---- Simulation parameter vector (trial-wise) ----
%     sim_params_vec = [ ...
%         sim.dt, ...
%         fixed.sigma_sq, ...
%         fixed.x0, ...
%         theta_i, ...
%         lambda_i ];
% 
%     % ---- Unique RNG seed per trial ----
%     unique_seed   = uint64(idx_bout * 1000 + randi([1 60]));
%     trial_seed(idx_bout) = unique_seed;
% 
%     % ---- Run high-resolution DDM simulation ----
%     res = sim_ddm_seeded(fine_signal, sim_params_vec, 1, unique_seed);

    mu_ig = theta_s ./ mu_s;
    lambda_ig = theta_s .^ 2;

    % Generate one sample from inverse Gaussian
    rt.decision(idx_trials) = random('InverseGaussian', mu_ig, lambda_ig);
    rt.ig(idx_trials) = rt.decision(idx_trials) + tndt_s;

end

toc

rt = [rt ncomp_vars];

% Plot the two resulting distributions, they should match

fh = figure('color', 'w', 'Position', [100, 100, 600, 600]);
tiledlayout(2, 1, 'TileSpacing', 'compact', 'Padding', 'loose')

nexttile
hold on
histogram(rt.st, 0:1/30:sim_params.T, 'EdgeColor', 'none', 'Normalization', 'pdf')
histogram(rt.ig, 0:1/30:sim_params.T, 'EdgeColor', 'none', 'Normalization', 'pdf')
apply_generic(gca)
xlabel('Duration (s)'); ylabel('pdf')
xlim([0 2])
ylim([0 1.5])

nexttile
histogram(abs(rt.st - rt.ig), 0:1/10:sim_params.T, 'EdgeColor', 'none', 'Normalization', 'pdf', 'FaceColor', 'k')
apply_generic(gca)
xlabel('|t_{simul} - t_{sampled}| (s)'); ylabel('pdf')  

% For some simulation it will go past ddm_params.T, so we need to censor those.

%  Now we added our vector column to the bouts table.yo
bouts_proc.durations_s = rt.ig;

tic 
model_results = run_fitting_updated(bouts_proc, points, '1', paths, 'export', true, 'extra', [], 'ground_truth', gt_table, 'n_bads', 8, 'bads_display', 'none');

model_results = run_fitting_updated(bouts_proc, points, '2', paths, 'export', true, 'extra', [], 'ground_truth', gt_table, 'n_bads', 8, 'bads_display', 'none');
model_results = run_fitting_newer(bouts_proc, points, 'dddm2', paths, 'export', true, 'extra', [], 'ground_truth', gt_table, 'n_bads', 8, 'bads_display', false);
toc

%%
[fh, ax, ax_inset] = fd_conditions('results', model_results, 'no_y', true, 'vis', 'on');
overlay_fits(fh, ax, ax_inset, 'results', model_results, 'export', true, 'extra', extra);
%plot_fit('freezes', bouts_proc, 'results', model_results)
%plot_fit('results', model_results, 'conditions', true)

plot_estimates('results', model_results)

%%
model_results = run_fitting_newer(bouts_proc, points, 'dddm3', paths, 'export', true, 'extra', [], 'ground_truth', gt_table);
%plot_fit('freezes', bouts_proc, 'results', model_results)
%plot_fit('results', model_results, 'conditions', true)

plot_estimates('results', model_results)
