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
    struct('name', 'sm')}}, ...
    'ground_truth', 60, ...
    'link', link_linear ...
    );

% For theta 1
model.theta = struct(...
    'predictors', {{ ...
    struct('name', 'fs'), ...
    struct('name', 'ls'), ...
    struct('name', 'intercept') ...
    }}, ...
    'ground_truth', [0.3 -0.2 2], ...
    'link', link_linear ...
    );

% % Non decision time
% model.tndt = struct( ...
%     'predictors', {{ ...
%     struct('name', 'intercept') ...
%     }}, ...
%     'ground_truth', 0.2, ...
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

tic
for idx_trials = 1:sim_params.n_trials

    ons = bouts_proc.onsets(idx_trials);

    sum_motion = motion_cache(bouts_proc.fly(idx_trials));
    sm_raw{idx_trials} = sum_motion(ons:ons + chunk_len - 1);
    sm_chunk = sm_raw{idx_trials};

    mu_tv = gt_table.mu_sm .* sm_chunk;% + gt_table.mu1_intercept;
    mu_st = ncomp_vars.mu(idx_trials);
    theta_s = ncomp_vars.theta(idx_trials);
%    tndt_s = ncomp_vars.tndt(idx_trials);

    % Simulate RT from full DDM
    [rt.ed(idx_trials), traj_ed] = extrema_detection_new('mu_t', mu_tv, 'theta', theta_s, ...
        'z', sim_params.z, 'dt', sim_params.dt, 'T', sim_params.T);%, 'ndt', tndt_s);

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
model_results = run_fitting_newer(bouts_proc, points, 'ed4', paths, 'export', true, 'extra', extra, 'ground_truth', gt_table);
plot_estimates('results', model_results, 'export', true, 'paths', paths)

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