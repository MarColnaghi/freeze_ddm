% Single DDM test

clear all

% Load the table first. We will take advantage of an already existing
% dataset.
threshold_imm = 3; threshold_mob = 3; threshold_pc = 4; id_code = sprintf('imm%d_mob%d_pc%d', threshold_imm, threshold_mob, threshold_pc);
paths = path_generator('folder', 'fitting_tests/sddm', 'bouts_id', id_code);
load(fullfile(paths.dataset, 'bouts.mat'));
bouts_proc = data_parser_new(bouts, 'type', 'immobility', 'period', 'loom', 'window', 'le');

% Only save useful variables in the table
y = table;
y.sm = bouts_proc.avg_sm_freeze_norm;
y.sm = bouts_proc.avg_sm_freeze_norm;
y.fs = bouts_proc.avg_fs_1s_norm;
y.ln = bouts_proc.nloom_norm;
y.ls = bouts_proc.sloom_norm;
y.intercept = ones(height(y),1);   % N‑by‑1 column of ones
predictors = y.Properties.VariableNames;

ncomp_vars = table();
link_linear = @(x) x;     % log link for bound height
link_logistic = @(x) 1./(1 + exp(-x));     % log link for bound height

% For mu
model.mu = struct( ...
    'predictors', {{ ...
    struct('name', 'sm'), ...
    struct('name', 'intercept') ...
    }}, ...
    'mask', [1 0 0 0 1], ...
    'ground_truth', [0.9 0.2], ...
    'link', link_linear ...
    );

% For theta
model.theta = struct( ...
    'predictors', {{ ...
    struct('name', 'fs'), ...
    struct('name', 'ls'), ...
    struct('name', 'intercept') ...
    }}, ...
    'mask', [0 1 0 1 1], ...
    'ground_truth', [0.4 0.4 1.2], ...
    'link', link_linear ...
    );

% Non decision time
model.tndt = struct( ...
    'predictors', {{ ...
    struct('name', 'intercept') ...
    }}, ...
    'mask', [0 0 0 0 1], ...
    'ground_truth', 0.25, ...
    'link', link_linear ...
    );

x = [0.9 0.2 0.4 0.4 0.2 0.25];
ncomp_vars = evaluate_model(model, x, y);
[gt, lbl] = get_ground_truth_vector(model);
gt_table = array2table(gt, 'VariableNames', lbl);

% Specify the seed
rng(1);

% General simulation parameters
sim_params.n_trials = height(bouts_proc);
sim_params.dt = 1/600;
sim_params.T = 60;
sim_params.time_vector = 0:sim_params.dt:sim_params.T;
sim_params.z = 0;

% Simulation settings
sim_params.kde_grid = 0:1/600:120;
sim_params.eval_trials = sim_params.n_trials;
sim_params.num_sims = 20000;

% Censoring/Truncation
points.truncation = [];
points.censoring = sim_params.T;

% Initialize outputs
rt = table;
rt.st = nan(sim_params.n_trials, 1);
rt.ig = nan(sim_params.n_trials, 1);

% Extract Social Motion TimeSeries
chunk_len = length(sim_params.time_vector);

tic
for idx_trials = 1:sim_params.n_trials

    mu_s = ncomp_vars.mu(idx_trials);
    theta_s = ncomp_vars.theta(idx_trials);
    tndt_s = ncomp_vars.tndt(idx_trials);

    % Simulate RT from full DDM
    [rt.st(idx_trials), traj_st] = drift_diff_new('mu', mu_s, 'theta', theta_s, ...
        'z', sim_params.z, 'dt', sim_params.dt, 'T', sim_params.T, 'ndt', tndt_s);


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
histogram(rt.st, 1/120:1/60:sim_params.T, 'EdgeColor', 'none', 'Normalization', 'pdf')
histogram(rt.ig, 1/120:1/60:sim_params.T, 'EdgeColor', 'none', 'Normalization', 'pdf')

apply_generic(gca)
xlabel('Duration (s)'); ylabel('pdf')
xlim([0 2])
ylim([0 1])

nexttile
histogram(abs(rt.st - rt.ig), 0:.25:sim_params.T, 'EdgeColor', 'none', 'Normalization', 'pdf', 'FaceColor', 'k')
apply_generic(gca)
xlabel('|t_{simul} - t_{sampled}| (s)'); ylabel('pdf')

% For some simulation it will go past ddm_params.T, so we need to censor those.

rt.st(isnan(rt.st)) = sim_params.T + 1; 
points.censoring = sim_params.T;
points.truncation = 0.4;

%  Now we added our vector column to the bouts table.
bouts_proc.durations_s = rt.st;

model_results = run_fitting_newer(bouts_proc, points, 'sddm1', paths, 'export', true, 'extra', [], 'ground_truth', gt_table);
plot_fit('results', model_results)
plot_estimates('results', model_results)

model_results = run_fitting_newer(bouts_proc, points, 'sddm2', paths, 'export', true, 'extra', [], 'ground_truth', gt_table);
plot_fit('results', model_results)
plot_estimates('results', model_results)

% run_fitting_newer(bouts_proc, points, 'sddm2', paths, 'export', true, 'extra', [], 'ground_truth', gt_table);
% run_fitting_newer(bouts_proc, points, 'sddm3', paths, 'export', true, 'extra', [], 'ground_truth', gt_table);
% run_fitting_newer(bouts_proc, points, 'sddm4', paths, 'export', true, 'extra', [], 'ground_truth', gt_table);
% run_fitting_newer(bouts_proc, points, 'sddm5', paths, 'export', true, 'extra', [], 'ground_truth', gt_table);