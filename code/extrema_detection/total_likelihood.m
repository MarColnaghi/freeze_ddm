clearvars

true_drift_scale = 1;
true_lambda = 0;
true_bound = 2;

% -- Likelihood / Signal Parameters (Coarse) --
fixed.dt = 1/60;       % 1ms for signal and PDE
fixed.dx = 0.01;
fixed.sigma_sq = 1.0;
%     fixed.bound = 1.0;
fixed.x0 = 0.0;

% -- Simulation Parameters (Fine) --
sim.dt = 0.00005;       % 0.05ms for ground truth simulation

% PDE Grid Setup
x_min = -7;
grid_size = round((true_bound - x_min) / fixed.dx) + 1;
start_idx = round((fixed.x0 - x_min) / fixed.dx) + 1;

% Save for objective function
fixed.x_min = x_min;
fixed.grid_size = grid_size;
fixed.start_idx = start_idx;

% Truncation and T_max
fixed.T_trunc = 0.3;
fixed.T_max = 10.5 - 1/60;

% Load the table first. We will take advantage of an already existing
% dataset.
threshold_imm = 2; threshold_mob = 2; threshold_pc = 4; id_code = sprintf('imm%d_mob%d_pc%d', threshold_imm, threshold_mob, threshold_pc);
model_str = 'dddm2';
model_func = str2func(strcat('model_',model_str));
model_struct = model_func();

col = cmapper();

color_list = {col.timevarying_sm, col.extremadetection};
model_list = {'dddm2', 'ded2'};
run_list = {'run13', 'run10'};
label_list = {'Integration', 'Extrema Detection'};

paths = path_generator('folder', 'fitting_freezes/le');

results_folder = fullfile(paths.results, model_list{1}, run_list{1});
temp_results = importdata(fullfile(results_folder, sprintf('fit_results_%s.mat', model_list{1})));
surrogate = importdata(fullfile(results_folder, 'surrogate.mat'));

for idx_model = 1:length(model_list)

    results_folder = fullfile(paths.results, model_list{idx_model}, run_list{idx_model});

    model_results = importdata(fullfile(results_folder, sprintf('fit_results_%s.mat', model_list{idx_model})));
    freezes = importdata(fullfile(results_folder, 'surrogate.mat'));
    extra = importdata(fullfile(results_folder, 'extra.mat'));

    est_params = table2array(model_results.estimates_mean(:, ~ismissing(model_results.estimates_mean)));
    ec = extra;
    ec.tndt = 0;
    ec = rmfield(ec, 'tndt');

    for idx_bout = 1:height(freezes)
        
        freeze = freezes(idx_bout, :);
        ec.soc_mot_array = extra.soc_mot_array(idx_bout, :);
        soc_mot_cell = {ec.soc_mot_array};

        ll(idx_bout, idx_model) =  - nll_fly_ddm_newer(est_params, freeze, model_results.points, strcat('model_', model_results.fitted_model), 'iid', '', ec);

        l_ma(idx_bout) = nll_fly_ddm_newer([1 2 0], freeze, model_results.points, strcat('model_', 'ed0'), 'iid', '', ec);
        l_ju(idx_bout) = log_joint_likelihood_3D([1 1e-12 2], soc_mot_cell, freeze.durations_s, fixed);

        [l_ma, l_ju]
    end
end


[sorted_ll, sorted_idx] = sort(diff(ll, [], 2));

%%


for idx_model = 1:length(model_list)
    i = 0;
    hold on

    results_folder = fullfile(paths.results, model_list{idx_model}, run_list{idx_model});

    model_results = importdata(fullfile(results_folder, sprintf('fit_results_%s.mat', model_list{idx_model})));
    freezes = importdata(fullfile(results_folder, 'surrogate.mat'));
    extra = importdata(fullfile(results_folder, 'extra.mat'));
    freezes.durations_s(freezes.durations_s > model_results.points.censoring) = 11;

    est_params = table2array(model_results.estimates_mean(:, ~ismissing(model_results.estimates_mean)));
    ec = extra;
    ec.tndt = 0;
    ec = rmfield(ec, 'tndt');

    for idx_bout = sorted_idx(end -4:end)'

        i = i + 1;
        fh = figure(i);
        fh.Position = [100 100 750 550];
        fh.Color = 'w';

        freeze = freezes(idx_bout, :);
        ec.soc_mot_array = extra.soc_mot_array(idx_bout, :);

        soc_mot_cell = {ec.soc_mot_array};
        [~, f, fd] =  nll_fly_ddm_newer(est_params, freeze, model_results.points, strcat('model_', model_results.fitted_model), 'iid', 'p', ec);

        ll = log_joint_likelihood_3D([1 2 0], soc_mot_cell, freeze.durations_s, fixed);

        plot(fd, f, 'LineWidth', 1.5, 'Color',  color_list{idx_model}, 'LineStyle', '--')
        plot(model_results.points.censoring, f(end), 'o', 'LineWidth', 1, 'MarkerSize', 5, 'MarkerEdgeColor', color_list{idx_model});

        xline(freeze.durations_s, 'k--', 'LineWidth', 2)
        apply_generic(gca, 'ylim', [0 2])

        %plot(fd(2:end), ec.soc_mot_array)
    end
end


function ll = log_joint_likelihood_3D(params, drifts_cell, rts, fixed)
prop_mean_drift = params(1);
prop_lambda     = params(2);
prop_bound      = params(3); % NEW

% --- Dynamic Grid Setup ---
% We must calculate the grid parameters for THIS specific bound.
% x_min is fixed, dx is fixed.

% Calculate number of spatial nodes needed to reach this bound
prop_grid_size = round((prop_bound - fixed.x_min) / fixed.dx) + 1;

% Calculate where x0 (0.0) falls on this grid
prop_start_idx = round((fixed.x0 - fixed.x_min) / fixed.dx) + 1;

% Pack MEX params
mex_p = [fixed.dt, fixed.dx, fixed.sigma_sq, fixed.x0, ...
    fixed.x_min, prop_grid_size, prop_start_idx, prop_lambda];

N = length(rts);
log_liks = zeros(N, 1);

t_vec = (0:fixed.dt:fixed.T_max)';

idx_trunc = round(fixed.T_trunc / fixed.dt);

parfor i = 1:N
    current_drift = drifts_cell{i} * prop_mean_drift;

    % 1. Solve PDE (Get PDF and Final Survival Prob)
    [p_dist, survival_at_Tmax] = leaky_pde_robust(current_drift, mex_p);

    % 2. Calculate Normalization Factor (The Truncation Logic)
    % We need the total probability mass that crossed BEFORE T_trunc.
    % This is the area under the curve from 0 to T_trunc.
    if idx_trunc > 1
        % Sum flux up to truncation point
        prob_missing = sum(p_dist(1:idx_trunc)) * fixed.dt;
    else
        prob_missing = 0;
    end

    % The total probability of a valid observation in the truncated window
    normalization_factor = 1.0 - prob_missing;

    % Safety: If the model predicts EVERYTHING happens before T_trunc,
    % the normalization approaches 0.
    if normalization_factor < 1e-12
        normalization_factor = 1e-12;
    end

    % 3. Handle Observed Data

    % A. Censored Data (Subject didn't respond by Tmax)
    if rts(i) >= (fixed.T_max - fixed.dt) || isnan(rts(i))
        % Likelihood is Survival Probability renormalized
        lik = survival_at_Tmax / normalization_factor;

        % B. Valid Response Data
    else
        idx = round(rts(i) / fixed.dt);

        % Sanity check: If data violates truncation (impossible data)
        if idx < idx_trunc
            % This data point exists but is below the truncation threshold
            % defined in the likelihood. This implies a mismatch definition.
            lik = 1e-12;
        else
            % Standard lookup
            if idx > length(p_dist), idx = length(p_dist); end

            %                 raw_pdf = p_dist(idx);
   
            raw_pdf = interp1(t_vec,p_dist,rts(i));

            % Renormalize PDF
            lik = raw_pdf / normalization_factor;
        end
    end

    % 4. Numerical Stability
    if lik < 1e-12 || isnan(lik)
        lik = 1e-12;
    end

    log_liks(i) = log(lik);
end

ll = sum(log_liks);
end

