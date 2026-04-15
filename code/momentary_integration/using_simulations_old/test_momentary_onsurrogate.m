clear all;

for idx_seed = 15
    
    exporting = true;
    col = cmapper();
    sim_params.rng = idx_seed;
    rng(sim_params.rng);

    % Model
    model = 'dddm2';
    select_run = 'run03';

    % Load important files
    paths = path_generator('folder', fullfile('/fitting_freezes/le', model, select_run));
    load(fullfile(paths.results, 'surrogate.mat'))
    load(fullfile(paths.results, sprintf('fit_results_%s.mat', model)));
    gt = table2array(estimates_mean(:, find(~ismissing(estimates_mean))));
    bouts_proc = surrogate;

    % Load the motion ts
    sim_params.motion_cache_path = fullfile(paths.dataset, 'motion_cache.mat');
    load(sim_params.motion_cache_path)

    % Immidiately switch paths for saving outputs accordingly
    paths = path_generator('folder', fullfile('/momentary_integration/surrogate', model));
    mkdir(paths.results); mkdir(paths.fig);
    create_output_dirs(paths)

    % General DDM parameters
    sim_params.dt = 1/60;
    sim_params.T = points.censoring - gt(end);
    sim_params.time_vector = 0:sim_params.dt:sim_params.T;
    sim_params.z = 0;

    % Simulation settings
    sim_params.kde_grid = 0:1/60:30;
    sim_params.eval_trials = height(bouts_proc);
    sim_params.num_sims = 25000;
    sim_params.norm_fact = 10;
    sim_params.n_trials = height(bouts_proc);

    % Create table
    y = table;
    y.sm = bouts_proc.avg_sm_freeze_norm;
    y.fs = bouts_proc.avg_fs_1s_norm;
    y.ln = bouts_proc.nloom_norm;
    y.ls = bouts_proc.sloom_norm;
    y.intercept = ones(height(y),1);
    y.smp = bouts_proc.avg_sm_freeze_norm;
    y.onsets = bouts_proc.onsets;
    y.fly = bouts_proc.fly;
    y.id = bouts_proc.id;

    predictors = y.Properties.VariableNames;
    model_eval = eval(sprintf('model_%s', model));
    ncomp_vars = evaluate_model(model_eval, gt, y);

    % Initialize outputs
    rt = table;
    rt_tv = nan(sim_params.n_trials, 1);
    rt_st = nan(sim_params.n_trials, 1);
    rt_ig = nan(sim_params.n_trials, 1);
    trial_type = nan(sim_params.n_trials, 1);

    sm_avg = nan(sim_params.n_trials, 1);
    sm_raw = cell(sim_params.n_trials, 1);
    timelock_sm = cell(sim_params.n_trials, 1);
    all_traj_tv = cell(sim_params.n_trials, 1);
    all_traj_st = cell(sim_params.n_trials, 1);
    all_traj_tv_det = cell(sim_params.n_trials, 1);
    all_traj_tv_sto = cell(sim_params.n_trials, 1);
    % Extract Social Motion TimeSeries
    chunk_len = length(sim_params.time_vector);

    tic
    for idx_trials = 1:sim_params.n_trials

        ons = bouts_proc.onsets(idx_trials);

        sum_motion = motion_cache(bouts_proc.fly(idx_trials));
        sm_raw{idx_trials} = sum_motion(ons:ons + chunk_len - 1);
        sm_chunk = sm_raw{idx_trials};

        % Determine model 1 or 2 based on pmix

        if rand <= ncomp_vars.pmix(idx_trials)
            trial_type(idx_trials) = 1;
            mu_t = estimates_mean.mu1_sm .* sm_chunk;
            theta_s = ncomp_vars.theta1(idx_trials);
            mu_s = ncomp_vars.mu1(idx_trials);

        else
            trial_type(idx_trials) = 2;
            mu_t = estimates_mean.mu2_sm .* sm_chunk;
            theta_s = ncomp_vars.theta2(idx_trials);
            mu_s = ncomp_vars.mu2(idx_trials);

        end

        tndt_s = ncomp_vars.tndt(idx_trials);

        % Time-varying DDM
        [rt_tv(idx_trials), traj_tv, t, traj_tv_det, traj_tv_sto] = drift_diff_new('mu_t', mu_t, 'theta', theta_s, ...
            'z', sim_params.z, 'dt', sim_params.dt, 'T', sim_params.T, 'ndt', tndt_s);

        % Stationary DDM based on average SM until response
        steps_tv = round((rt_tv(idx_trials) - tndt_s) ./ sim_params.dt);
        sm_avg(idx_trials) = mean(sm_chunk(1:min(chunk_len, steps_tv)));

        if trial_type(idx_trials) == 1
            mu_avg = estimates_mean.mu1_sm *  sm_avg(idx_trials);
        else
            mu_avg = estimates_mean.mu2_sm *  sm_avg(idx_trials);
        end

        % Stationary DDM based on average SM until response
        [rt_st(idx_trials), traj_st] = drift_diff_new('mu', mu_avg, 'theta', theta_s, ...
            'z', sim_params.z, 'dt', sim_params.dt, 'T', sim_params.T, 'ndt',  tndt_s);

        % Store results
        timelock_sm{idx_trials} = sm_chunk(1:min(chunk_len, steps_tv));

        all_traj_tv{idx_trials} = traj_tv;
        all_traj_tv_sto{idx_trials} = traj_tv_sto;
        all_traj_tv_det{idx_trials} = traj_tv_det;
        all_traj_st{idx_trials} = traj_st;

        mu_ig = theta_s ./ mu_avg;
        lambda_ig = theta_s .^ 2;

        % Generate one sample from inverse Gaussian
        rt_decision = random('InverseGaussian', mu_ig, lambda_ig);
        rt_ig(idx_trials) = rt_decision + tndt_s;

    end
    toc

    % This Takes around 7s.

    % Record results in the table
    rt.tv = rt_tv; rt.st = rt_st; rt.ig = rt_ig; rt.sm_avg = sm_avg; rt.sm_ts = sm_raw; rt.id = bouts_proc.id;
    diff_rt = rt.tv - rt.st;
    max(rt.tv)
    max(rt.st)

    % Let's make sure to save the results
    cd(paths.results)
    save('rt.mat', 'rt')
    save('sim_params.mat', 'sim_params')

    % Here we do some plots.
    endpoint = points.censoring + 1/20;
    fh = figure('color','w','Position',[100, 100, 800, 800]);
    tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact')
    nexttile
    hold on
    histogram(rt.tv, 1/120:1/20:endpoint, 'Normalization', 'pdf', 'EdgeColor', 'none', 'FaceColor', col.timevarying_sm, 'FaceAlpha', 0.6);
    histogram(rt.st, 1/120:1/20:endpoint, 'Normalization', 'pdf', 'EdgeColor', 'none', 'FaceColor', col.stationary_sm, 'FaceAlpha', 0.6);
    xlabel('Freeze Time (s)'); ylabel('pdf');
    ax = gca;
    apply_generic(ax, 24);
    xlim([0 11])
    ylim([-0.03 1.53])
    xline(points.censoring, 'LineWidth', .5, 'LineStyle', '-.', 'Label', sprintf('rt > t_{cens}: tv: %d - st: %d', sum(isnan(rt.tv)), sum(isnan(rt.st))), 'LabelOrientation', 'aligned',  'FontSize', 14)

    nexttile
    histogram(diff_rt, -10:1/20:10, 'Normalization', 'pdf', 'EdgeColor', 'none', 'FaceColor', 'k', 'FaceAlpha', 0.6);
    xline(mean(diff_rt, 'omitnan'), 'r--')
    xlabel('Delta Freeze Time (s)'); ylabel('pdf');
    ax = gca;
    apply_generic(ax, 24);

    nexttile
    hold on
    histogram(rt.ig, 1/120:1/20:endpoint, 'Normalization', 'pdf', 'EdgeColor', 'none', 'FaceColor', col.closedform, 'FaceAlpha', 0.6);
    histogram(rt.st, 1/120:1/20:endpoint, 'Normalization', 'pdf', 'EdgeColor', 'none', 'FaceColor', col.stationary_sm, 'FaceAlpha', 0.6);
    xlabel('Freeze Time (s)'); ylabel('pdf');
    ax = gca;
    apply_generic(ax, 24);
    xlim([0 11])
    ylim([-0.03 1.53])
    xline(points.censoring, 'LineWidth', .5, 'LineStyle', '-.', 'Label', sprintf('rt > t_{cens}: ig: %d', sum(rt.ig > points.censoring)), 'LabelOrientation', 'aligned',  'FontSize', 14)

    y.durations_s = rt.ig;
    [~, f, fd] = nll_fly_ddm_newer(gt, y, points, sprintf('model_%s', model), 'iid', 'p', []);
    nexttile
    hold on
    histogram(rt.ig(rt.ig >= 0.5), 1/120 + tndt_s:1/20:max(rt.ig), 'Normalization', 'pdf', 'EdgeColor', 'none', 'FaceColor', col.closedform, 'FaceAlpha', 0.6);
    plot(fd, f, 'k--')
    xlabel('Freeze Time (s)'); ylabel('pdf');
    ax = gca;
    apply_generic(ax, 24);
    xlim([0 11])
    exporter(fh, paths, 'rts_summary.pdf', 'flag', exporting)

    % Plot offset-locked sm
    [sorted_len, sort_idx] = sort(cellfun(@length, timelock_sm), 'descend');
    aligned_mat = padder_for_imagesc(timelock_sm(sort_idx), sorted_len, 'offset');
    fh = plot_sta_sm(aligned_mat, 'direction', 'offset', 'clim', [0 1], 'center', 'mean');
    exporter(fh, paths, 'offset_aligned_rts.pdf', 'flag', exporting)

    %% Here we fit the models to the two datasets

    points.truncation = [];
    rt.st(isnan(rt.st)) = 11;
    rt.tv(isnan(rt.tv)) = 11;
    rt.ig(rt.ig > points.censoring) = 11;

    % Here is important to modify the table to add sm only to the 'sm' column.
    % pmix should not be affected.

    y.sm = rt.sm_avg;

    % Which datasets you would like to test
    tested_datasets = {'tv', 'st', 'ig'};

    %%%%%%% Loop FIT %%%%%%%

    for gen_type = tested_datasets

        % Select dataset and modify durations column in the table
        gen_data = gen_type{1};
        y.durations_s = rt.(gen_data);

        % No need to create path, we will have a global mat file with all the fits.
        paths_temp.fig = fullfile(paths.fig, sprintf('fit_%s', gen_data));
        mkdir(paths_temp.fig)

        % Fit the dataset
        model_results.(gen_data) = run_fitting_newer(y, points, model, paths, 'export', false, 'extra', [], 'ground_truth', estimates_mean);

        % Just present two simple fits
        plot_fit('freezes', y, 'results', model_results.(gen_data), 'export', exporting, 'paths', paths_temp)
        plot_estimates('results', model_results.(gen_data), 'export', exporting, 'paths', paths_temp)


    end

    cd(paths.results)
    save('model_results.mat', 'model_results')

    model_results.(gen_data).points
    %%%%%%% Loop LIKELIHOOD %%%%%%%

    for gen_type = tested_datasets

        % Select dataset and modify durations column in the table
        gen_data = gen_type{1};
        y.durations_s = rt.(gen_data);

        paths_temp.results = fullfile(paths.results, sprintf('sims_%s', gen_data));
        mkdir(paths_temp.results); cd(paths_temp.results)

        est_params = table2array(model_results.(gen_data).estimates_mean(:, find(~ismissing(model_results.(gen_data).estimates_mean))));
        lls_output = evaluate_likelihood_parallel(y, motion_cache, sim_params, model_results.(gen_data), []);

        save('lls_final.mat', 'lls_output')
        save('y.mat', 'y')
        save(sprintf('all_traj_%s.mat', gen_data), sprintf('all_traj_%s', gen_data))

    end

    
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
