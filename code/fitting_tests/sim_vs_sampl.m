% Double DDM test

clear all

% Load the table first. We will take advantage of an already existing
% dataset.
threshold_imm = 0; threshold_mob = 0; threshold_pc = 4; id_code = sprintf('imm%d_mob%d_pc%d', threshold_imm, threshold_mob, threshold_pc);
paths = path_generator('folder', 'fitting_tests/dddm', 'bouts_id', id_code);
load(fullfile(paths.dataset, 'bouts.mat'));
bouts_proc = data_parser_new(bouts, 'type', 'immobility', 'period', 'loom', 'window', 'le');

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
model.mu1 = struct( ...
    'predictors', {{ ...
    struct('name', 'intercept')}}, ...
    'ground_truth', 1.3, ...
    'link', link_linear ...
    );

% For theta 1
model.theta1 = struct( ...
    'predictors', {{struct('name', 'intercept')
    }}, ...
    'ground_truth', 4, ...
    'link', link_linear ...
    );

% Non decision time
model.tndt = struct( ...
    'predictors', {{ ...
    struct('name', 'intercept') ...
    }}, ...
    'ground_truth', 0.2, ...
    'link', link_linear ...
    );

[gt, lbl] = get_ground_truth_vector(model);
x = gt(~isnan(gt));
gt_table = array2table(gt, 'VariableNames', lbl);
ncomp_vars = evaluate_model(model, gt_table, y);

% Specify the seed
sim_params.rng = 15;
rng(sim_params.rng);

% General simulation parameters
sim_params.n_trials = height(bouts_proc);
sim_params.dt = 1/600;
sim_params.T = 300;
sim_params.time_vector = 0:sim_params.dt:sim_params.T;
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
rt.st = nan(sim_params.n_trials, 1);
rt.ig = nan(sim_params.n_trials, 1);
trial_type = nan(sim_params.n_trials, 1);

% Extract Social Motion TimeSeries
chunk_len = length(sim_params.time_vector);

tic
for idx_trials = 1:sim_params.n_trials

    mu_s = ncomp_vars.mu1(idx_trials);
    theta_s = ncomp_vars.theta1(idx_trials);
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
histogram(rt.st, 0:1/5:sim_params.T, 'EdgeColor', 'none', 'Normalization', 'pdf')
histogram(rt.ig, 0:1/5:sim_params.T, 'EdgeColor', 'none', 'Normalization', 'pdf')
xlim([0 30])
apply_generic(gca)