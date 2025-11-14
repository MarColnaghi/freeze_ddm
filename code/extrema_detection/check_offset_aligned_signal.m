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
    struct('name', 'sm') ...
    }}, ...
    'ground_truth', 2.2, ...
    'link', link_linear ...
    );

% For theta 1
model.theta = struct(...
    'predictors', {{ ...
    struct('name', 'intercept') ...
    }}, ...
    'ground_truth', 0.39, ...
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
for idx_trials = 1:sim_params.n_trials

    ons = bouts_proc.onsets(idx_trials);

    sum_motion = motion_cache(bouts_proc.fly(idx_trials)) ./ 10;
    sm_raw{idx_trials} = sum_motion(ons:ons + chunk_len - 1);
    sm_chunk = sm_raw{idx_trials};

    mu_tv = gt_table.mu_sm .* sm_chunk;
    mu_st = ncomp_vars.mu(idx_trials);
    theta_s = ncomp_vars.theta(idx_trials);
    tndt_s = ncomp_vars.tndt(idx_trials);

    % Simulate RT from full DDM
    [rt.ed(idx_trials), traj_ed] = extrema_detection_new('mu_t', mu_tv, 'theta', theta_s, ...
        'z', sim_params.z, 'dt', sim_params.dt, 'T', sim_params.T, 'ndt', samples_sec(idx_trials)); % samples_sec(idx_trials));

    t_of_cross = round(rt.ed(idx_trials) * 1/sim_params.dt);

    t_of_cross = min(chunk_len, t_of_cross);
    timelock_sm{idx_trials} = sum_motion(ons:ons + t_of_cross - 1);
end

toc
timelock_sm{isnan(rt.ed)} = {};
rt.ed(rt.ed > points.censoring) = nan;
rt.ed(isnan(rt.ed)) = nan;
points.censoring = sim_params.T;
points.truncation = min(rt.ed) - 1/60;

fh = figure('color', 'w', 'Position', [100, 100, 600, 300]);
[sorted_len, sort_idx] = sort(cellfun(@length, timelock_sm), 'descend');
aligned_mat = padder_for_imagesc(timelock_sm(sort_idx), sorted_len, 'offset');
fh = plot_sta_sm(aligned_mat, 'direction', 'offset', 'clim', [0 1], 'center', 'mean');