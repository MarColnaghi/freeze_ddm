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
model = model_ed3();

% Specify Ground Truth
model.mu.ground_truth = 10;
model.theta.ground_truth = [0.3 0.3];
model.tndt.ground_truth = 0;

% Specify the seed
sim_params.rng = 9999;
rng(sim_params.rng);

% General simulation parameters
sim_params.n_trials = height(bouts_proc);
sim_params.dt = 1/60;
sim_params.T = 10.5;
sim_params.time_vector = 0:sim_params.dt:sim_params.T;
sim_params.z = 0;

fps = 60;
sim_params.exp_time_vector = 0:1/fps:sim_params.T;
points.censoring = sim_params.T;
points.truncation = 0;

% Initialize outputs
rt = table;
rt.ed = nan(sim_params.n_trials, 1);
trial_type = nan(sim_params.n_trials, 1);

sm_during = extract_sm_from_bouts(bouts_proc, 'type', 'onsets', 'output_type', 'mat', 'window', [0 630]);

sm_during_us = interp1(sim_params.exp_time_vector, sm_during', sim_params.time_vector, 'nearest')';
bouts_proc.sm = sm_during_us;

[gt, lbl] = get_ground_truth_vector(model);
x = gt(~isnan(gt));
gt_table = array2table(gt, 'VariableNames', lbl);
ncomp_vars = evaluate_model(model, gt_table, bouts_proc);

% Extract Social Motion TimeSeries
chunk_len = length(sim_params.time_vector);

tic
for idx_trials = 1:height(bouts_proc)

    mu_tv = ncomp_vars.mu(idx_trials, :);
    %mu_st = ncomp_vars.mu(idx_trials);
    theta_s = ncomp_vars.theta(idx_trials);
    tndt_s = ncomp_vars.tndt(idx_trials);

    % Simulate RT from full DDM
    [rt.ed(idx_trials), traj_ed] = extrema_detection_new('mu_t', mu_tv, 'theta', theta_s, ...
        'z', sim_params.z, 'dt', sim_params.dt, 'T', sim_params.T, 'ndt', tndt_s); % samples_sec(idx_trials));

end
toc

rt.ed(isnan(rt.ed)) = 20; 
points.censoring = sim_params.T;
points.truncation = 0;

fh = figure('color', 'w', 'Position', [100, 100, 600, 300]);
tiledlayout(1, 1, 'TileSpacing', 'compact', 'Padding', 'loose')
nexttile
hold on
histogram(rt.ed, 1/120:1/20:sim_params.T + sim_params.dt * 3, 'EdgeColor', 'none', 'Normalization', 'pdf')
apply_generic(gca)
ylabel('Density'); xlabel('Duration (s)'); ax.YAxis.Visible = 'off';
exporter(fh, paths, 'Durations.pdf')

%%
extra.soc_mot_array = sm_during_us;

model_2_fit = 'ed3';
model_results = run_fitting_newer(bouts_proc, points, model_2_fit, paths, 'export', true,  'ground_truth', gt_table, 'bads_display', 'iter', 'pass_ndt', true, 'n_bads', 3, 'extra', extra);
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