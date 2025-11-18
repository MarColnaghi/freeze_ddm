clear all

% Load the table first. We will take advantage of an already existing
% dataset.
threshold_imm = 3; threshold_mob = 3; threshold_pc = 4; id_code = sprintf('imm%d_mob%d_pc%d', threshold_imm, threshold_mob, threshold_pc);
paths = path_generator('folder', 'extrema_detection/visual_aid', 'bouts_id', id_code);
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
    struct('name', 'sm')}}, ...
    'ground_truth', [0.7], ...
    'link', link_linear ...
    );

% For theta 1
model.theta = struct(...
    'predictors', {{ ...
    struct('name', 'intercept') ...
    }}, ...
    'ground_truth', [0.80], ...
    'link', link_linear ...
    );

% Non decision time
model.tndt = struct( ...
    'predictors', {{ ...
    struct('name', 'intercept') ...
    }}, ...
    'ground_truth', 0, ...
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

for idx_trials = 333%:sim_params.n_trials

    ons = bouts_proc.onsets(idx_trials);

    sum_motion = motion_cache(bouts_proc.fly(idx_trials));
    sm_raw{idx_trials} = sum_motion(ons:ons + chunk_len - 1) ./ 10;
    sm_chunk = sm_raw{idx_trials};

    mu_tv = gt_table.mu_sm .* sm_chunk;
    mu_st = ncomp_vars.mu(idx_trials);
    theta_s = ncomp_vars.theta(idx_trials);
    tndt_s = ncomp_vars.tndt(idx_trials);

    % Simulate RT from full DDM
    [rt.ed(idx_trials), traj_ed] = extrema_detection_new('mu_t', mu_tv * 30, 'theta', theta_s, ...
        'z', sim_params.z, 'dt', sim_params.dt, 'T', sim_params.T, 'ndt', tndt_s); % samples_sec(idx_trials));

    [rt.st(idx_trials), traj_st] = drift_diff_new('mu_t', mu_tv, 'theta', theta_s, ...
        'z', sim_params.z, 'dt', sim_params.dt, 'T', sim_params.T, 'ndt', tndt_s); % samples_sec(idx_trials));

end


fh = figure('color', 'w', 'Position', [100, 100, 400, 300]);
tiledlayout(1, 1, 'TileSpacing', 'compact', 'Padding', 'loose')
nexttile
hold on
apply_generic(gca)
plot(0:length(traj_st) - 1, traj_st, 'b-', 'LineWidth', 1.4)

imagesc(1:length(mu_tv), -1.5, mu_tv')
colormap(cbrewer2('Reds'))

xlim([0 70]);
ylim([-1.3 1.3])
yline(x(2), 'k-',  'LineWidth', 2.5)
yticks([])
xlabel('samples');
exporter(fh, paths, 'int.pdf')

fh = figure('color', 'w', 'Position', [100, 100, 400, 300]);
tiledlayout(1, 1, 'TileSpacing', 'compact', 'Padding', 'loose')
nexttile
hold on
apply_generic(gca)
plot(0:length(traj_ed) - 1, traj_ed, 'ro', 'LineWidth', 1.4)
imagesc(1:length(mu_tv), -1.5, mu_tv')
colormap(cbrewer2('Reds'))
xlim([0 70]);
ylim([-1.3 1.3])
yline(x(2), 'k-',  'LineWidth', 2.5)
yticks([])
xlabel('samples');
exporter(fh, paths, 'ed.pdf')
