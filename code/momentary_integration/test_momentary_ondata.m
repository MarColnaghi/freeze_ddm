clear all
idx_seed = randi(100);

exporting = true;
col = cmapper();
sim_params.rng = idx_seed;
rng(sim_params.rng);

% Model
model = 'dddm2';
select_run = 'run03';
gen_data = 'fr';

% Load important files
paths = path_generator('folder', fullfile('/fitting_freezes/le', model, select_run));
load(fullfile(paths.results, 'surrogate.mat'))
model_results.(gen_data) = load(fullfile(paths.results, sprintf('fit_results_%s.mat', model)));
est_params = table2array(model_results.(gen_data).estimates_mean(:, find(~ismissing(model_results.(gen_data).estimates_mean))));
bouts_proc = surrogate;

% Load the motion ts
sim_params.motion_cache_path = fullfile(paths.dataset, 'motion_cache.mat');
load(sim_params.motion_cache_path)

% Immidiately switch paths for saving outputs accordingly
paths = path_generator('folder', fullfile('/momentary_integration/freezes', model));
mkdir(paths.results); mkdir(paths.fig);
create_output_dirs(paths)

% General DDM parameters
sim_params.dt = 1/60;
sim_params.T = model_results.(gen_data).points.censoring - est_params(end);
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
y.durations_s = bouts_proc.durations_s;

predictors = y.Properties.VariableNames;
model_eval = eval(sprintf('model_%s', model));
ncomp_vars = evaluate_model(model_eval, est_params, y);

y.durations_s(y.durations_s > model_results.(gen_data).points.censoring, :) = 11;

rt = table(); rt.(gen_data) = y.durations_s; rt.id = y.id;

save('rt.mat', 'rt') 
save('model_results.mat', 'model_results');
save('sim_params.mat', 'sim_params') 

paths_temp.results = fullfile(paths.results, sprintf('sims_%s', gen_data));
mkdir(paths_temp.results); cd(paths_temp.results)
lls_output = evaluate_likelihood_parallel(y, motion_cache, sim_params, model_results.(gen_data), []);

save('lls_final.mat', 'lls_output')
save('y.mat', 'y')


% % Store in output
% lls.(gen_data).lls_timevarying             = lls_output(:,1);
% lls.(gen_data).lls_stationary              = lls_output(:,2);
% lls.(gen_data).lls_stationary_cf           = lls_output(:,3);
% 
% [idx.(gen_data)] = plots_lls(lls.(gen_data));
% 
% select_indexes = idx.(gen_data)([1, 2, 3, end - 2, end - 1, end]);
% lls_output = evaluate_likelihood(synth_data, rt, est_params, ddm_params, model, points, col, select_indexes);



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
    cd(paths.results)
    
    % Also update figure path to match the new run
    paths.fig = fullfile(paths.fig, run_name);
    mkdir(paths.fig);

    % Assign the updated paths back to base workspace (if needed)
    assignin('caller', 'paths', paths);
end
