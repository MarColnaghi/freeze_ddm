close all; clear all;
%rng(1);

doplot = true;
col = cmapper();
lls = table;

% Model
model = 'dddm2';
select_run = 'run02';

% Here you should load the estimated parameters.
paths = path_generator('folder', fullfile('/fitting_freezes/le', model, select_run));
load(fullfile(paths.results, 'surrogate.mat'))
load(fullfile(paths.results, sprintf('fit_results_%s.mat', model)));
est_params = table2array(estimates_mean(:, find(~ismissing(estimates_mean))));
bouts_proc = surrogate;

% Simulation settings
sim_params.kde_grid = 0:1/60:30;
sim_params.eval_trials = height(bouts_proc);
sim_params.num_sims = 25000;
sim_params.norm_fact = 10;

% General DDM parameters
sim_params.dt = 1/300;
sim_params.T = 59;
sim_params.time_vector = 0:sim_params.dt:sim_params.T;
sim_params.z = 0;

bouts_proc = bouts_proc(bouts_proc.onsets < 36000 - sim_params.T * 60, :);
sim_params.n_trials = height(bouts_proc);

load(fullfile(paths.dataset, 'motion_cache.mat'))

y = table;
y.sm = bouts_proc.avg_sm_freeze_norm;
y.sm = bouts_proc.avg_sm_freeze_norm;
y.fs = bouts_proc.avg_fs_1s_norm;
y.ln = bouts_proc.nloom_norm;
y.ls = bouts_proc.sloom_norm;
y.intercept = ones(height(y),1);   % N‑by‑1 column of ones
predictors = y.Properties.VariableNames;
model_eval = eval(sprintf('model_%s', model));
ncomp_vars = evaluate_model(model_eval, est_params, y);

% Initialize outputs
rt = table;
rt_tv = nan(sim_params.n_trials, 1);
rt_st = nan(sim_params.n_trials, 1);
rt_ed = nan(sim_params.n_trials, 1);
rt_ig = nan(sim_params.n_trials, 1);
trial_type = nan(sim_params.n_trials, 1);

sm_avg = nan(sim_params.n_trials, 1);

timelock_sm = cell(sim_params.n_trials, 1);
all_traj = cell(sim_params.n_trials, 1);

% Extract Social Motion TimeSeries
chunk_len = length(sim_params.time_vector);

tic
for idx_trials = 1:sim_params.n_trials

    ons = bouts_proc.onsets(1);

    sum_motion = motion_cache(bouts_proc.fly(idx_trials));
    sm_raw{idx_trials} = sum_motion(ons:ons + chunk_len - 1);
    
    sm_chunk = sm_raw{idx_trials};

    % Determine model 1 or 2 based on pmix
    if rand < ncomp_vars.pmix(idx_trials)
        
        trial_type(idx_trials) = 1;
        mu_t = ncomp_vars.mu1(idx_trials) .* sm_chunk;
        theta_s = ncomp_vars.theta1(idx_trials);
        
    else
        trial_type(idx_trials) = 2;
        mu_t = ncomp_vars.mu2(idx_trials) .* sm_chunk;
        theta_s = ncomp_vars.theta2(idx_trials);

    end
    
    tndt_s = ncomp_vars.tndt(idx_trials);

    % Time-varying DDM
    [rt_tv(idx_trials), traj_tv] = drift_diff_new('mu_t', mu_t, 'theta', theta_s, ...
        'z', sim_params.z, 'dt', sim_params.dt, 'T', sim_params.T, 'ndt', tndt_s);

    %[rt_ed(idx_trials), traj_ed] = extrema_detection_new('mu_t', mu_t, 'theta', theta_s, ...
   %     'dt', sim_params.dt, 'T', sim_params.T, 'ndt', sim_params.ndt);

    if isnan(rt_tv(idx_trials))
       disp('over censoring')
    end

    % Stationary DDM based on average SM until response
    steps_tv = round((rt_tv(idx_trials) - tndt_s) ./ sim_params.dt);
    sm_avg(idx_trials) = mean(sm_chunk(1:min(chunk_len, steps_tv)));
    
    if trial_type(idx_trials) == 1 
        mu_avg = estimates_mean.mu1_sm *  sm_avg(idx_trials);
    else
        mu_avg = estimates_mean.mu2_sm *  sm_avg(idx_trials);
    end

    % Stationary DDM based on average SM until response
    %steps_ed = round((rt_ed(idx_trials) - sim_params.ndt) * 60);
    
    [rt_st(idx_trials), traj_st] = drift_diff_new('mu', mu_avg, 'theta', theta_s, ...
        'z', sim_params.z, 'dt', sim_params.dt, 'T', sim_params.T, 'ndt',  tndt_s);

    % Store results
    timelock_sm{idx_trials} = sm_chunk(1:min(chunk_len, steps_tv));
    %timelock_ed{idx_trials} = sm_chunk(1:max(1, steps_ed));

    all_traj{idx_trials} = traj_tv;

    mu_ig = theta_s ./ mu_avg;
    lambda_ig = theta_s .^ 2;

    % Generate one sample from inverse Gaussian
    rt_decision = random('InverseGaussian', mu_ig, lambda_ig);
    rt_ig(idx_trials) = rt_decision + tndt_s;

end
toc

% This Takes around 7s.

% Record results in the table
rt.tv = rt_tv; rt.st = rt_st; rt.ig = rt_ig; rt.sm_avg = sm_avg; rt.sm_ts = sm_raw';
diff_rt = rt.tv - rt.st;

fh = figure('color','w','Position',[100, 100,900, 400]);
tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact')
nexttile
hold on
histogram(rt.tv, 1/120 + tndt_s:1/20:sim_params.T, 'Normalization', 'pdf', 'EdgeColor', 'none', 'FaceColor', col.timevarying_sm, 'FaceAlpha', 0.6);
histogram(rt.st, 1/120 + tndt_s:1/20:sim_params.T, 'Normalization', 'pdf', 'EdgeColor', 'none', 'FaceColor', col.stationary_sm, 'FaceAlpha', 0.6);
xlabel('Freeze Time (s)'); ylabel('pdf');
ax = gca;
apply_generic(ax, 24);
xlim([0 11])

nexttile
hold on
histogram(rt.ig, 1/120 + tndt_s:1/20:sim_params.T, 'Normalization', 'pdf', 'EdgeColor', 'none', 'FaceColor', col.timevarying_sm, 'FaceAlpha', 0.6);
histogram(rt.st, 1/120 + tndt_s:1/20:sim_params.T, 'Normalization', 'pdf', 'EdgeColor', 'none', 'FaceColor', col.stationary_sm, 'FaceAlpha', 0.6);
xlabel('Freeze Time (s)'); ylabel('pdf');
ax = gca;
apply_generic(ax, 24);
xlim([0 11])