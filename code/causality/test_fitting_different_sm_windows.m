% Double DDM test

clearvars

% Load the table first. We will take advantage of an already existing
% dataset.
threshold_imm = 2; threshold_mob = 2; threshold_pc = 4; id_code = sprintf('imm%d_mob%d_pc%d', threshold_imm, threshold_mob, threshold_pc);
paths = path_generator('folder', 'causality/fitting_windows', 'bouts_id', id_code, 'imfirst', false);
bouts = importdata(fullfile(paths.dataset, 'bouts.mat'));
motion_cache = importdata(fullfile(paths.cache_path, 'motion_cache.mat'));

bouts_proc = data_parser_new(bouts, 'type', 'immobility', 'period', 'loom', 'window', 'le');
points.censoring = 10.5;
points.truncation = 0.3;

%  Now we added our vector column to the bouts table

chunk_len = points.censoring * 60;
sm_ts = nan(height(bouts_proc), chunk_len);
sm_raw = cell(1 - height(bouts_proc));

for idx_trials = 1:height(bouts_proc)

    ons = bouts_proc.onsets(idx_trials);
    off = bouts_proc.ends(idx_trials) - 1;
    sum_motion = motion_cache(bouts_proc.fly(idx_trials));
    sm_raw{idx_trials} = sum_motion(ons:ons + chunk_len - 1) ./ 10;

    sm_ts(idx_trials, :) = sum_motion(ons - 630:ons - 1) ./ 10;

end

model = 'dddm2';
paths = path_generator('folder', sprintf('causality/fitting_windows/%s', model), 'bouts_id', id_code, 'imfirst', false);
create_output_dirs(paths);

for idx_backframes = [600]

    backframes = idx_backframes;

    mkdir(fullfile(paths.fig, num2str(backframes)));
    mkdir(fullfile(paths.results, num2str(backframes)));
    paths_backframe.results = fullfile(paths.results, num2str(backframes));
    paths_backframe.fig = fullfile(paths.fig, num2str(backframes));

    bouts_proc.avg_sm_pre_norm = mean(sm_ts(:, end - backframes:end), 2);
    soc_mot_array = cell2mat(sm_raw)';
    extra.soc_mot_array = soc_mot_array;
    model_results = run_fitting_newer(bouts_proc, points, model, paths, 'export', false,...
        'bads_display', false, 'pass_ndt', false, 'n_bads', 2, 'extra', extra);


    model_results.bouts_path = paths_backframe.results; model_results.paths = paths.fig;

    model_results.motion_cache_path = fullfile(paths.cache_path, 'motion_cache.mat');

    save(fullfile(paths_backframe.results, sprintf('fit_results_%s.mat', model)), '-struct', 'model_results');
    save(fullfile(paths_backframe.results, 'extra.mat'), 'extra');
    save(fullfile(paths_backframe.results, 'surrogate.mat'), 'bouts_proc');

end

% %%
% fh = plot_fit('results', model_results, 'conditions', false, 'export', true, 'bin_size', 3, 'censored_inset', true, 'type', 'continuous');
% fh_conditions = plot_fit('results', model_results, 'conditions', true, 'export', true, 'bin_size', 10 , 'type', 'continuous');
% plot_estimates('results', model_results, 'export', true, 'paths', paths, 'ylimits', [-1 4])
% 


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