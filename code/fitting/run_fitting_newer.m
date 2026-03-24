function [model_results, estimates] = run_fitting_newer(freeze_table, points, idx_model, paths, varargin)

% RUN_FITTING_NEWER - Run DDM model fitting with BADS and VBMC
% Version: Refactored July 2025

%% Parse Optional Inputs
opt = inputParser;
addParameter(opt, 'extra', []);
addParameter(opt, 'export', false);
addParameter(opt, 'bads_display', 'iter');
addParameter(opt, 'pass_ndt', false);
addParameter(opt, 'iid', false);
addParameter(opt, 'ground_truth', []);
addParameter(opt, 'only_bads', false);
addParameter(opt, 'n_bads', 5);
addParameter(opt, 'vbmc_exhaustive', false);

parse(opt, varargin{:});

extra = opt.Results.extra;
export_results = opt.Results.export;
ground_truth = opt.Results.ground_truth;
bads_display = opt.Results.bads_display;
pass_ndt = opt.Results.pass_ndt;
n_bads = opt.Results.n_bads;
vbmc_exhaustive = opt.Results.vbmc_exhaustive;
only_bads = opt.Results.only_bads;

%% Truncation Filter

if isfield(points, 'truncation') && ~isempty(points.truncation)
    mask = freeze_table.durations_s >= points.truncation;
    freeze_table = freeze_table(mask, :);
    if isfield(extra, 'soc_mot_array')
        extra.soc_mot_array = extra.soc_mot_array(mask, :);
    end
    fprintf('The truncation point is %.2fs \n', points.truncation);
else
    fprintf('There is no truncation point \n');

end

fprintf('The censoring point is %.2fs \n', points.censoring);
fprintf('The smallest rt is %.2fs \n', min(freeze_table.durations_s));

%% Model Setup
% Extract the Bounds
model_str = sprintf('model_%s', idx_model);
model_func = str2func(model_str);
model_obj = model_func();

% Extract the Bounds
[LB, PLB, PUB, UB, ~] = extract_bounds_from_model(model_obj);
[~, lbl, mask] = get_ground_truth_vector(model_obj);

%% Setting up the table

if ismember('smp', freeze_table.Properties.VariableNames)
    freeze_table.smp = freeze_table.smp;
else
    freeze_table.smp = freeze_table.sm;
end

freeze_table.intercept = ones(height(freeze_table), 1);

%% BADS Optimization

% Likelihood Function
llfun = @(x) nll_fly_ddm_newer(x, freeze_table, points, model_str, 'iid', 'n', extra);

num_iters = n_bads;
options_bads = bads('defaults');
options_bads.Display = bads_display;

nvars = numel(PLB);
eval_param = zeros(num_iters, nvars);
fval = zeros(num_iters, 1);

tic
all_estimates = nan(num_iters, length(lbl));
res_table = array2table(all_estimates, 'VariableNames', lbl);

fprintf('\n--- Starting BADS Multi-Run Optimization ---\n');

for idx = 1:num_iters
    x0 = PLB + rand(1, nvars) .* (PUB - PLB);

    [eval_param(idx,:), fval(idx)] = bads(llfun, x0, LB, UB, PLB, PUB, [], options_bads);

    current_full_row = nan(1, length(lbl));
    current_full_row(mask) = eval_param(idx, :);

    res_table{idx, :} = current_full_row;
    is_empty_col = all(isnan(res_table{idx, :}), 1);

    clean_table = res_table(:, ~is_empty_col);
    disp(clean_table);

end

is_empty_col = all(isnan(res_table{:, :}), 1);
clean_table = res_table(:, ~is_empty_col);

meta_data = table((1:num_iters)', fval, 'VariableNames', {'Iter', 'NLL'});
clean_table = [meta_data, clean_table];

disp('--- Final Optimization Results ---');
disp(clean_table);
toc

[~, best_idx] = sort(fval);
eval_param = eval_param(best_idx,:); fval = fval(best_idx, :);

estimates = nan(1, length(lbl));
estimates(find(mask)) = eval_param(1, :);
temp_table = array2table(estimates, 'VariableNames', lbl);

if only_bads
    estimates = eval_param(1, :);
    model_results = [];
else


    if pass_ndt
        extra.tndt = temp_table.tndt_intercept;
        eval_param = eval_param(:, 1:end-1);
        LB = LB(1:end-1);
        PLB = PLB(1:end-1);
        PUB = PUB(1:end-1);
        UB = UB(1:end-1);

    end


    % Find fixed parameters
    is_fixed = (LB == PLB) & (PLB == PUB) & (PUB == UB);

    % Optionally store fixed values
    extra.fixed_param = struct();

    idx_params = find(mask);

    if any(is_fixed)
        fixed_names  = lbl(idx_params(is_fixed));
        fixed_values = LB(is_fixed);

        for i = 1:numel(fixed_names)
            extra.fixed_param.(fixed_names{i}) = fixed_values(i);
        end
    end

    % Remove fixed parameters from optimization vectors
    eval_param(:, is_fixed) = [];
    LB(is_fixed)  = [];
    PLB(is_fixed) = [];
    PUB(is_fixed) = [];
    UB(is_fixed)  = [];

    if ~isempty(ground_truth)
        if width(ground_truth) == width(array2table(estimates, 'VariableNames', lbl))
            ground_truth = [ground_truth; array2table(estimates, 'VariableNames', lbl)];
        else
            try
                ground_truth = outerjoin(ground_truth, array2table(estimates, 'VariableNames', lbl), 'MergeKeys', true);
            catch
                disp('not sharing any parameter')
                temp = array2table(estimates, 'VariableNames', lbl);
                temp.(ground_truth.Properties.VariableNames{1}) = nan;
                ground_truth = outerjoin(ground_truth, temp, 'MergeKeys', true);
            end
        end
        plc_hold = ground_truth;
        plc_hold(:, all(ismissing(plc_hold))) = [];
        disp(plc_hold)
    else
        starting_point = array2table(estimates, 'VariableNames', lbl);
        starting_point(:, ismissing(starting_point)) = [];
        disp(starting_point)
    end

    %% VBMC Optimization
    llfun = @(x) -nll_fly_ddm_newer(x, freeze_table, points, model_str, 'iid', 'n', extra);
    lpriorfun = @(x) structured_prior(x, prior_info, LB, PLB, PUB, UB);

    lpriorfun = @(x) msplinetrapezlogpdf(x, LB, PLB, PUB, UB);
    postfun = @(x) lpostfun(x, llfun, lpriorfun);

    options_vbmc.Display = 'iter';
    options_vbmc.MaxFunEvals = 800;

    if vbmc_exhaustive
        options_vbmc.SpecifyTargetNoise = false;
        options_vbmc.TolStableCount = 80;
        options_vbmc.MinFinalComponents = 50;

        % 1. Force a high-resolution initial map
        options_vbmc.FunEvalStart = 5;
        % 2. Force it to run longer (don't stop on stability early)
        options_vbmc.MinFunEvals = 250;
        % 3. Force higher precision before declaring convergence
        % options_vbmc.TolImprovement = 0.004;
        % 4. Start with a more complex mixture model (smoother tails)
        options_vbmc.Kwarmup = 5;
    end

    [VP, ELBO, ELBO_SD] = vbmc(postfun, eval_param(1,:), LB, UB, PLB, PUB, options_vbmc);
    [x_mean, x_sigma] = vbmc_moments(VP);
    x_std = sqrt(diag(x_sigma));
    vbmc_plot(VP);

    %% Store Fit Results
    model_results = struct;
    model_results.elbo = [ELBO, ELBO_SD];
    model_results.elbo_normalized = [ELBO, ELBO_SD] ./ height(freeze_table);
    model_results.time = datetime;

    [~, lbl, mask] = get_ground_truth_vector(model_obj);
    estimates_mean = nan(1, length(lbl));
    estimates_std = nan(1, length(lbl));

    if pass_ndt
        x_mean = [x_mean extra.tndt];
        x_std  = [x_std; 0];
    end

    estimates_mean(find(mask)) = x_mean;
    estimates_std(find(mask)) = x_std;

    model_results.estimates_mean = array2table(estimates_mean, 'VariableNames', lbl);
    model_results.estimates_std = array2table(estimates_std, 'VariableNames', lbl);

    model_results.points.truncation = points.truncation;
    model_results.points.censoring = points.censoring;

    model_results.starting_position = eval_param(1,:);

    model_results.fitted_model = idx_model;
    model_results.vp = VP;

    if ~isempty(ground_truth)
        model_results.ground_truth = ground_truth;

    end

    %plot_params_new(model, paths, model_result)

    %% Export Results
    if export_results

        %% Prepare Paths
        paths.fig = fullfile(paths.fig, idx_model);
        paths.results = fullfile(paths.results, idx_model);
        create_output_dirs(paths);
        model_results.bouts_path = paths.results; model_results.fig_path = paths.fig;

        model_results.motion_cache_path = fullfile(paths.cache_path, 'motion_cache.mat');

        save(fullfile(paths.results, 'model_results.mat'), '-struct', 'model_results');
        save(fullfile(paths.results, 'freeze.mat'), 'freeze_table');
        save(fullfile(paths.results, 'extra.mat'), 'extra');

    end
end
end

function create_output_dirs(paths)
    if ~exist(paths.fig, 'dir'), mkdir(paths.fig); end
    if ~exist(paths.results, 'dir'), mkdir(paths.results); end

    % Use a temporary variable to look for existing runs in the BASE directory
    run_folders = dir(fullfile(paths.results, 'run*'));
    run_nums = [];
    for i = 1:length(run_folders)
        if run_folders(i).isdir
            tokens = regexp(run_folders(i).name, '^run(\d+)', 'tokens');
            if ~isempty(tokens)
                run_nums(end+1) = str2double(tokens{1}{1});
            end
        end
    end

    if isempty(run_nums)
        next_run = 1;
    else
        next_run = max(run_nums) + 1;
    end

    datestamp = string(datetime('now', 'Format', 'yyMMdd'));
    run_name = sprintf('run%02d_%s', next_run, datestamp);

    % Create the specific run subfolders
    new_results_path = fullfile(paths.results, run_name);
    new_fig_path = fullfile(paths.fig, run_name);

    mkdir(new_results_path);
    mkdir(new_fig_path);

    % Update the paths struct for the caller
    paths.results = new_results_path;
    paths.fig = new_fig_path;
    
    assignin('caller', 'paths', paths);
end


function lp = structured_prior(x, prior_spec, LB, PLB, PUB, UB)

persistent call_count
if isempty(call_count)
    call_count = 0;
end
call_count = call_count + 1;

plot_every = 50;   % <-- adjust as needed

lp = 0;

for i = 1:numel(x)

    switch prior_spec(i).type

        case 'trapezoidal'
            lp = lp + msplinetrapezlogpdf( ...
                x(i), LB(i), PLB(i), PUB(i), UB(i));

        case 'exponential'
            if x(i) < LB(i) || x(i) > UB(i)
                lp = -Inf;
                return
            end
            lambda = prior_spec(i).lambda;
            lp = lp + log(lambda) - lambda * x(i);

        otherwise
            error('Unknown prior type: %s', prior_spec(i).type)
    end
end

% ---- Diagnostic plotting (occasionally)
if mod(call_count, plot_every) == 0
    plot_prior_snapshot(prior_spec, LB, PLB, PUB, UB, x);
end

end

function plot_prior_snapshot(prior_spec, LB, PLB, PUB, UB, x_curr)

npar = numel(prior_spec);
clf
tiledlayout('flow')

for i = 1:npar
    nexttile; hold on

    xgrid = linspace(LB(i), UB(i), 300);
    logp = nan(size(xgrid));

    switch prior_spec(i).type
        case 'exponential'
            lambda = prior_spec(i).lambda;
            logp = log(lambda) - lambda * xgrid;

        case 'trapezoidal'
            logp = arrayfun(@(x) ...
                msplinetrapezlogpdf(x, LB(i), PLB(i), PUB(i), UB(i)), ...
                xgrid);
    end

    % Normalize for visualization
    p = exp(logp - max(logp));
    p = p / trapz(xgrid, p);

    plot(xgrid, p, 'k', 'LineWidth', 1.5)

    % --- Plausible bounds
    yl = ylim;
    plot([PLB(i) PLB(i)], yl, '--', 'Color',[0.7 0.7 0.7])
    plot([PUB(i) PUB(i)], yl, '--', 'Color',[0.7 0.7 0.7])

    % --- Current parameter value
    xc = x_curr(i);
    plot([xc xc], yl, 'r-', 'LineWidth', 2)

    % --- Visual warning if outside plausible region
    if xc < PLB(i) || xc > PUB(i)
        title(sprintf('param %d ⚠', i))
    else
        title(sprintf('param %d', i))
    end

    box on
end

drawnow limitrate
end



