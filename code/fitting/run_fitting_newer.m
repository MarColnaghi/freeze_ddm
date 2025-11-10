function [model_results] = run_fitting_newer(surrogate, points, idx_model, paths, varargin)

% RUN_FITTING_NEWER - Run DDM model fitting with BADS and VBMC
% Version: Refactored July 2025

%% Parse Optional Inputs
opt = inputParser;
addParameter(opt, 'extra', []);
addParameter(opt, 'export', false);
addParameter(opt, 'bads_display', false);
addParameter(opt, 'pass_ndt', true);

addParameter(opt, 'ground_truth', []);
addParameter(opt, 'only_bads', false);

parse(opt, varargin{:});

extra = opt.Results.extra;
export_results = opt.Results.export;
ground_truth = opt.Results.ground_truth;
bads_display = opt.Results.bads_display;
pass_ndt = opt.Results.pass_ndt;

%% Truncation Filter
if isfield(points, 'truncation') && ~isempty(points.truncation)
    mask = surrogate.durations_s >= points.truncation;
    surrogate = surrogate(mask, :);
    if isfield(extra, 'soc_mot_array')
        extra.soc_mot_array = extra.soc_mot_array(mask, :);
    end
    fprintf('The truncation point is %.2fs \n', points.truncation);
else
    fprintf('There is no truncation point \n');

end

fprintf('The censoring point is %.2fs \n', points.censoring);
fprintf('The smallest rt is %.2fs \n', min(surrogate.durations_s));
fprintf('We are trying to recover a model \n');

%% Model Setup
model_str = sprintf('model_%s', idx_model);
model_out = eval(model_str);
[LB, PLB, PUB, UB] = extract_bounds_from_model(model_out);

%% Objective Function
surrogate.sm = surrogate.avg_sm_freeze_norm;
surrogate.smp = surrogate.avg_sm_freeze_norm;
surrogate.fs = surrogate.avg_fs_1s_norm;
surrogate.ln = surrogate.nloom_norm;
surrogate.ls = surrogate.sloom_norm;
surrogate.intercept = ones(height(surrogate), 1);
llfun = @(x) nll_fly_ddm_newer(x, surrogate, points, model_str, 'iid', 'n', extra);

%% BADS Optimization
num_iters = 2;
if bads_display
    options_bads.Display = 'iter';
else
    options_bads.Display = 'none';
end
nvars = numel(PLB);
x0_all = PLB + rand(num_iters, nvars) .* (PUB - PLB);

eval_param = zeros(num_iters, nvars);
fval = zeros(num_iters, 1);
tic
for idx = 1:num_iters
    fprintf('Currently bads run #%d \n', idx)
    [eval_param(idx,:), fval(idx)] = bads(llfun, x0_all(idx,:), LB, UB, PLB, PUB, [], options_bads);
    eval_param(idx,:)
end

toc
[~, best_idx] = sort(fval);
eval_param = eval_param(best_idx,:); fval = fval(best_idx, :);

[~, lbl, mask] = get_ground_truth_vector(model_out);
estimates = nan(1, length(lbl));
estimates(find(mask)) = eval_param(1, :);
temp_table = array2table(estimates, 'VariableNames', lbl);

if pass_ndt
    extra.tndt = temp_table.tndt_intercept;
    eval_param = eval_param(:, 1:end-1);
    LB = LB(1:end-1);
    PLB = PLB(1:end-1);
    PUB = PUB(1:end-1);
    UB = UB(1:end-1);

end

if ~isempty(ground_truth)
    if width(ground_truth) == width(array2table(estimates, 'VariableNames', lbl))
        ground_truth = [ground_truth; array2table(estimates, 'VariableNames', lbl)];
    else
        ground_truth = outerjoin(ground_truth, array2table(estimates, 'VariableNames', lbl), 'MergeKeys', true);
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
llfun = @(x) -nll_fly_ddm_newer(x, surrogate, points, model_str, 'iid', 'n', extra);
lpriorfun = @(x) msplinetrapezlogpdf(x, LB, PLB, PUB, UB);
postfun = @(x) lpostfun(x, llfun, lpriorfun);

options_vbmc.Display = 'iter';
options_vbmc.MaxFunEvals = 500;
options_vbmc.RetryMaxFunEvals = 3 * 500;

[VP, ELBO, ELBO_SD] = vbmc(postfun, eval_param(1,:), LB, UB, PLB, PUB, options_vbmc);
[x_mean, x_sigma] = vbmc_moments(VP);
x_std = sqrt(diag(x_sigma));
vbmc_plot(VP);

%% Store Fit Results
model_results = struct;
model_results.elbo = [ELBO, ELBO_SD];
model_results.elbo_normalized = [ELBO, ELBO_SD] ./ height(surrogate);
model_results.time = datetime;

[~, lbl, mask] = get_ground_truth_vector(model_out);
estimates_mean = nan(1, length(lbl));
estimates_std = nan(1, length(lbl));

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
    model_results.bouts_path = paths.results;
    model_results.motion_cache_path = fullfile(paths.dataset, 'motion_cache.mat');
    
    save(fullfile(paths.results, sprintf('fit_results_%s.mat', idx_model)), '-struct', 'model_results');
    save(fullfile(paths.results, 'surrogate.mat'), 'surrogate');
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
