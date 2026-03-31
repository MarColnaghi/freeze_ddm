% Double DDM test

clearvars

% Open CMAPPER
col = cmapper();

% Load the table first. We will take advantage of an already existing
% dataset. 
threshold_imm = 2; threshold_mob = 2; threshold_pc = 4; id_code = sprintf('imm%d_mob%d_pc%d', threshold_imm, threshold_mob, threshold_pc);
paths = path_generator('folder', fullfile('fitting_tests','dddm'), 'bouts_id', id_code);
load(fullfile(paths.dataset, 'bouts.mat'));
bouts_proc = data_parser_new(bouts, 'type', 'immobility', 'period', 'loom', 'window', 'le');
motion_cache = importdata(fullfile(paths.cache_path, 'motion_cache.mat'));

ncomp_vars = table();
link_linear = @(x) x;     % log link for bound height
link_logistic = @(x) 1./(1 + exp(-x));     % log link for bound height

% Pick a model
model = model_dddm2();

% Specify Ground Truth
model.mu1.ground_truth = 1.4;
model.theta1.ground_truth = 0.7;
model.mu2.ground_truth = 1.7;
model.theta2.ground_truth = [0 2.2];
model.pmix.ground_truth = [0 0 0 0.2];
model.tndt.ground_truth = 0.3;

% Specify the seed
sim_params.rng = 9999;
rng(sim_params.rng);

% General simulation parameters
sim_params.n_trials = height(bouts_proc);
sim_params.dt = 1/300;
sim_params.T = 10.5;
sim_params.time_vector = 0:sim_params.dt:sim_params.T;
sim_params.z = 0;

points.censoring = sim_params.T;
points.truncation = 0;

% Initialize outputs
rt = table;
rt.st = nan(sim_params.n_trials, 1);
rt.ig = nan(sim_params.n_trials, 1);
trial_type = nan(sim_params.n_trials, 1);

chunk_len = 630;
total_length = 30;

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

    mu_ig = theta_s ./ mu_s;
    lambda_ig = theta_s .^ 2;

    % Generate one sample from inverse Gaussian
    rt.decision(idx_trials) = random('InverseGaussian', mu_ig, lambda_ig);
    rt.ig(idx_trials) = rt.decision(idx_trials) + tndt_s;

end

toc
rt = [rt ncomp_vars];

% nan from the simulations are censored datapoints
rt.st(isnan(rt.st)) = 20;

% Plot the two resulting distributions, they should match

fh = figure('color', 'w', 'Position', [100, 100, 600, 600]);
tiledlayout(3, 1, 'TileSpacing', 'compact', 'Padding', 'loose')

nexttile
hold on
histogram(rt.st, 0:1/30:sim_params.T, 'EdgeColor', 'none', 'Normalization', 'pdf')
histogram(rt.ig, 0:1/30:sim_params.T, 'EdgeColor', 'none', 'Normalization', 'pdf')
apply_generic(gca)
xlabel('Duration (s)'); ylabel('pdf')
xlim([0 sim_params.T])
ylim([0 1.5])

nexttile
hold on
histogram(rt.st(trial_type == 1), 0:1/30:sim_params.T, 'EdgeColor', 'none', 'Normalization', 'pdf', 'FaceColor', col.processes.short)
histogram(rt.st((trial_type == 2)), 0:1/30:sim_params.T, 'EdgeColor', 'none', 'Normalization', 'pdf', 'FaceColor', col.processes.long)
apply_generic(gca)
xlabel('Duration (s)'); ylabel('pdf')
xlim([0 sim_params.T])
ylim([0 1.5])

nexttile
histogram(abs(rt.st - rt.ig), 0:1/10:sim_params.T, 'EdgeColor', 'none', 'Normalization', 'pdf', 'FaceColor', 'k')
apply_generic(gca)
xlabel('|t_{simul} - t_{sampled}| (s)'); ylabel('pdf')  

% For some simulation it will go past ddm_params.T, so we need to censor those.

%  Now we added our vector column to the bouts table. You can select which
%  of the two to add (ig vs st)

bouts_proc.durations_s = rt.st;

% And then fit the model
model_results = run_fitting_newer(bouts_proc, points, 'dddm15', paths, 'export', true, 'extra', [], 'ground_truth', gt_table, 'n_bads', 2, 'bads_display', 'none');

% Plot the posteriors distribution together with the ground truth
plot_estimates('results', model_results, 'export', true, 'ylimits', [-2 5])

[fh, ax, ax_inset] = fd_conditions('results', model_results, 'no_y', true, 'vis', true, 'bin_size', 3, 'col', 'gray2', 'censored_inset', true);
overlay_fits(fh, ax, ax_inset, 'results', model_results, 'export', true, 'extra', [], 'censored_inset', true, 'gt', true);
overlay_separate_processes(fh, ax, ax_inset, 'results', model_results, 'export', true, 'extra', [])
