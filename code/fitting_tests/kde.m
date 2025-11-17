clear all

% Load the table first. We will take advantage of an already existing
% dataset.
threshold_imm = 3; threshold_mob = 3; threshold_pc = 4; id_code = sprintf('imm%d_mob%d_pc%d', threshold_imm, threshold_mob, threshold_pc);
paths = path_generator('folder', 'fitting_tests/dddm', 'bouts_id', id_code);
load(fullfile(paths.dataset, 'bouts.mat'));
kde_estimates = importdata(fullfile('/Users/marcocolnaghi/PhD/freeze_ddm/model_results/fitting_freezes/bsl/kde_spontaneous', id_code, 'kde_estimates_bsl.mat'));
[~, idx] = unique(kde_estimates.Fkde, 'last');
Fkde = kde_estimates.Fkde(idx); xkde = kde_estimates.xkde(idx); fkde = kde_estimates.fkde(idx); 

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
model.mu = struct( ...
    'predictors', {{ ...
    struct('name', 'sm')}}, ...
    'ground_truth', 0.9, ...
    'link', link_linear ...
    );

% For theta 1
model.theta = struct( ...
    'predictors', {{struct('name', 'intercept')
    }}, ...
    'ground_truth', 1.4, ...
    'link', link_linear ...
    );

model.pmix = struct( ...
    'predictors', {{ ...
    struct('name', 'fs') ...
    struct('name', 'ln') ...
    }}, ...
    'ground_truth', [0.8 -0.4], ...
    'link', link_logistic ...
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
sim_params.rng = 165;
rng(sim_params.rng);

% General simulation parameters
sim_params.n_trials = height(bouts_proc);
sim_params.dt = 1/60;
sim_params.T = 10.5;
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
        mu_s = ncomp_vars.mu(idx_trials);
        theta_s = ncomp_vars.theta(idx_trials);
        tndt_s = ncomp_vars.tndt(idx_trials);

        % Simulate RT from full DDM
        [rt.st(idx_trials), traj_st] = drift_diff_new('mu', mu_s, 'theta', theta_s, ...
            'z', sim_params.z, 'dt', sim_params.dt, 'T', sim_params.T, 'ndt', tndt_s);

    else
        trial_type(idx_trials) = 2;
        U = rand;
        rt.st(idx_trials) = interp1(Fkde, xkde, U, 'linear', 'extrap');

    end

end
toc

rt = [rt ncomp_vars];
figure
histogram(rt.st(trial_type == 1), 1/120:1/20:300)
hold on
histogram(rt.st(trial_type == 2), 1/120:1/20:300)
xlim([0 10])

rt.st(isnan(rt.st)) = sim_params.T + 1; 
points.censoring = sim_params.T;
points.truncation = 0;

%  Now we added our vector column to the bouts table.
bouts_proc.durations_s = rt.st;
bouts_proc.sm = bouts_proc.avg_sm_freeze_norm;
bouts_proc.smp = bouts_proc.avg_sm_freeze_norm;
bouts_proc.fs = bouts_proc.avg_fs_1s_norm;
bouts_proc.ls = bouts_proc.sloom_norm;
bouts_proc.ln = bouts_proc.nloom_norm;
bouts_proc.intercept = ones(height(y),1);

extra.Fkde = Fkde;
extra.fkde = fkde;
extra.xkde = xkde;

model_results = run_fitting_newer(bouts_proc, points, 'ksddm1', paths, 'export', true, 'extra', extra, 'ground_truth', gt_table);

plot_estimates('results', model_results)
plot_fit('freezes', bouts_proc, 'results', model_results, 'extra', extra)
plot_fit('freezes', bouts_proc, 'results', model_results, 'extra', extra, 'conditions', true)
