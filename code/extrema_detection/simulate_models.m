clear all
close all

% Load the table first. We will take advantage of an already existing
% dataset.

col = cmapper();

threshold_imm = 3; threshold_mob = 3; threshold_pc = 4; id_code = sprintf('imm%d_mob%d_pc%d', threshold_imm, threshold_mob, threshold_pc);
paths = path_generator('folder', 'extrema_detection/ac_vs_ed', 'bouts_id', id_code);

load(fullfile(paths.dataset, 'bouts.mat'));
bouts_proc = data_parser_new(bouts, 'type', 'immobility', 'period', 'loom', 'window', 'le');
bouts_proc = bouts_proc(bouts_proc.nloom < 15, :);

% Load the motion ts
sim_params.motion_cache_path = fullfile(paths.dataset, 'motion_cache.mat');
load(sim_params.motion_cache_path)

% Only save useful variables in the table
y = table;
y.sm = bouts_proc.avg_sm_freeze_norm;
y.smp = bouts_proc.avg_sm_freeze_norm;
y.fs = bouts_proc.avg_fs_1s_norm;
y.ln = bouts_proc.nloom_norm;
y.ls = bouts_proc.sloom_norm;
y.intercept = ones(height(y),1);   % N‑by‑1 column of ones
predictors = y.Properties.VariableNames;

link_linear = @(x) x;     % log link for bound height
link_logistic = @(x) 1./(1 + exp(-x));     % log link for bound height

% For mu
model.mu = struct( ...
    'predictors', {{ ...
    struct('name', 'sm'), ...
    struct('name', 'intercept') ...
    }}, ...
    'mask', [1 0 0 0 0 1], ...
    'ground_truth', [1.2 0.4], ...
    'link', link_linear ...
    );

% For theta
model.theta = struct( ...
    'predictors', {{ ...
    struct('name', 'fs'), ...
    struct('name', 'ls'), ...
    struct('name', 'intercept') ...
    }}, ...
    'mask', [0 0 1 0 1 1], ...
    'ground_truth', [1 1 1], ...
    'link', link_linear ...
    );

% Non decision time
model.tndt = struct( ...
    'predictors', {{ struct('name', 'intercept') }}, ...
    'mask', [0 0 0 0 0 1], ...
    'ground_truth', 0, ...
    'link', link_linear ...
    );

[gt, lbl] = get_ground_truth_vector(model);
x = gt(~isnan(gt));
ncomp_vars = evaluate_model(model, x, y);
gt_table = array2table(gt, 'VariableNames', lbl);

% Specify the seed
rng(1);

% General simulation parameters
sim_params.n_trials = height(bouts_proc);
sim_params.dt = 1/60;
sim_params.T = 90;
sim_params.time_vector = sim_params.dt:sim_params.dt:sim_params.T;
sim_params.z = 0;
sim_params.snr = 100;
sim_params.gt_table = gt_table;

% Simulation settings
sim_params.kde_grid = 0:1/600:120;
sim_params.eval_trials = sim_params.n_trials;
sim_params.num_sims = 20000;

% Censoring/Truncation
sim_params.points.truncation = [];
sim_params.points.censoring = sim_params.T;

% Initialize outputs
rt = table;
rt_ac = nan(sim_params.n_trials, 1);
rt_ed = nan(sim_params.n_trials, 1);
sm_raw = cell(sim_params.n_trials, 1);

% Extract Social Motion TimeSeries
chunk_len = length(sim_params.time_vector);

tic
for idx_trials = 1:height(bouts_proc)
    
    ons = bouts_proc.onsets(idx_trials);

    sum_motion = motion_cache(bouts_proc.fly(idx_trials));
    sm_raw{idx_trials} = sum_motion(ons:ons + chunk_len - 1);
    sm_chunk = sm_raw{idx_trials};

    mu_tv = gt_table.mu_sm .* sm_chunk + gt_table.mu_intercept;
    theta_s = ncomp_vars.theta(idx_trials);
    tndt_s = ncomp_vars.tndt(idx_trials);

    % Simulate RT from full DDM
    [rt_ac(idx_trials), traj_ac] = drift_diff_new('mu_t', mu_tv, 'theta', theta_s, ...
        'z', sim_params.z, 'dt', sim_params.dt, 'T', sim_params.T, 'ndt', tndt_s, 'sigma', 1);

    [rt_ed(idx_trials), traj_ed] = extrema_detection_new('mu_t', mu_tv .* sim_params.snr, 'theta', theta_s, ...
        'z', sim_params.z, 'dt', sim_params.dt, 'T', sim_params.T, 'ndt', tndt_s, 'sigma', 1);

    
    if idx_trials < 50 & idx_trials > 40
        figure
        hold on
        plot(traj_ed)
        plot(traj_ac)
        plot(sm_chunk)
        hold on
    end
end

toc

rt.ac = rt_ac; rt.ed = rt_ed; rt.sm_ts = sm_raw;
fh = figure('color', 'w', 'Position', [100, 100, 600, 600]);

hold on
histogram(rt.ed, sim_params.dt/2:sim_params.dt * 5 :sim_params.T + 1, 'FaceColor', col.extremadetection, 'EdgeColor', 'none', 'Normalization', 'pdf')
histogram(rt.ac, sim_params.dt/2:sim_params.dt * 5 :sim_params.T + 1, 'FaceColor', col.timevarying_sm, 'EdgeColor', 'none', 'Normalization', 'pdf')
xline(mean(rt.ed, 'omitnan'), 'Color', col.extremadetection, 'LineWidth', 2)
xline(mean(rt.ac, 'omitnan'), 'Color', col.timevarying_sm, 'LineWidth', 2)
xlabel('Duration')
ylabel('Count')
apply_generic(gca, 24);

%
y = table;
y.sm = bouts_proc.avg_sm_freeze_norm;
y.smp = bouts_proc.avg_sm_freeze_norm;
y.fs = bouts_proc.avg_fs_1s_norm;
y.ln = bouts_proc.nloom_norm;
y.ls = bouts_proc.sloom_norm;
y.intercept = ones(height(y),1);
y.onsets = bouts_proc.onsets;
y.fly = bouts_proc.fly;
y.id = bouts_proc.id;

create_output_dirs(paths)

tested_datasets = {'ac', 'ed'};

%%%%%%% Loop FIT %%%%%%%

for gen_type = tested_datasets

    % Select dataset and modify durations column in the table
    gen_data = gen_type{1};
    y.durations_s = rt.(gen_data);

    y.durations_s(isnan(y.durations_s)) = sim_params.T + 1;
    % No need to create path, we will have a global mat file with all the fits.
    paths_temp.results = fullfile(paths.results, sprintf('sims_%s', gen_data));
    mkdir(paths_temp.results); cd(paths_temp.results)

    save('y.mat', 'y')

end


function create_output_dirs(paths)
    % Ensure base directories exist
    if ~exist(paths.fig, 'dir'), mkdir(paths.fig); end
    if ~exist(paths.results, 'dir'), mkdir(paths.results); end

    % Auto-incrementing run folder inside results
    run_folders = dir(fullfile(paths.results, 'run*'));
    run_nums = [];

    for i = 1:length(run_folders)
        if run_folders(i).isdir
            tokens = regexp(run_folders(i).name, '^run(\d+)$', 'tokens');
            if ~isempty(tokens)
                run_nums(end+1) = str2double(tokens{1}{1}); %#ok<AGROW>
            end
        end
    end

    if isempty(run_nums)
        next_run = 1;
    else
        next_run = max(run_nums) + 1;
    end

    run_name = sprintf('run%02d', next_run);
    paths.results = fullfile(paths.results, run_name);
    mkdir(paths.results);

    % Also update figure path to match the new run
    paths.fig = fullfile(paths.fig, run_name);
    mkdir(paths.fig);

    % Assign the updated paths back to base workspace (if needed)
    assignin('caller', 'paths', paths);
end



