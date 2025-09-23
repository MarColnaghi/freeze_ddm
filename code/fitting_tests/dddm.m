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
    struct('name', 'sm')}}, ...
    'mask', [1 0 0 0 0 0], ...
    'ground_truth', 0.9, ...
    'link', link_linear ...
    );

% For theta 1
model.theta1 = struct( ...
    'predictors', {{struct('name', 'intercept')
    }}, ...
    'mask', [0 0 0 0 1 0], ...
    'ground_truth', 0.8, ...
    'link', link_linear ...
    );

% For mu 2
model.mu2 = struct( ...
    'predictors', {{ ...
    struct('name', 'sm'), ...
    }}, ...
    'mask', [1 0 0 0 0 0], ...
    'ground_truth', [0.6], ...
    'link', link_linear ...
    );

% For theta
model.theta2 = struct( ...
    'predictors', {{ ...
    struct('name', 'ls'), ...
    struct('name', 'intercept') ...
    }}, ...
    'mask', [0 0 0 1 1 0], ...
    'ground_truth', [0.7 2.2], ...
    'link', link_linear ...
    );

model.pmix = struct( ...
    'predictors', {{ ...
    struct('name', 'sm'), ...
    struct('name', 'fs') ...
    struct('name', 'ln') ...
    struct('name', 'ls') ...
    }}, ...
    'mask', [1 1 1 1 0 0], ...
    'ground_truth', [2.0 -1.0 0.8 -0.4], ...
    'link', link_logistic ...
    );

% Non decision time
model.tndt = struct( ...
    'predictors', {{ ...
    struct('name', 'intercept') ...
    }}, ...
    'mask', [0 0 0 0 1 0], ...
    'ground_truth', 0.25, ...
    'link', link_linear ...
    );

x = [0.9 0.8 0.6 0.7 2.2 2.0 -1.0 0.8 -0.4 0.25];
ncomp_vars = evaluate_model(model, x, y);
[gt, lbl] = get_ground_truth_vector(model);
gt_table = array2table(gt, 'VariableNames', lbl);

% Specify the seed
rng(2);

% General simulation parameters
sim_params.n_trials = height(bouts_proc);
sim_params.dt = 1/1200;
sim_params.T = 60;
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

rt.st(isnan(rt.st)) = sim_params.T + 1; 
points.censoring = sim_params.T;
points.truncation = 0;

%  Now we added our vector column to the bouts table.
bouts_proc.durations_s = rt.st;
bouts_proc.sm = bouts_proc.avg_sm_freeze_norm;
bouts_proc.fs = bouts_proc.avg_fs_1s_norm;
bouts_proc.ls = bouts_proc.sloom_norm;
bouts_proc.ln = bouts_proc.nloom_norm;
bouts_proc.intercept = ones(height(y),1);
bouts_proc.smp = bouts_proc.avg_sm_freeze_norm;

model_results = run_fitting_newer(bouts_proc, points, 'dddm1', paths, 'export', true, 'extra', [], 'ground_truth', gt_table);
plot_fit('freezes', bouts_proc, 'results', model_results)
plot_fit('results', model_results, 'conditions', true)

plot_estimates('results', model_results)

model_results = run_fitting_newer(bouts_proc, points, 'dddm2', paths, 'export', true, 'extra', [], 'ground_truth', gt_table);
plot_fit('freezes', bouts_proc, 'results', model_results)
plot_fit('results', model_results, 'conditions', true)

plot_estimates('results', model_results)
