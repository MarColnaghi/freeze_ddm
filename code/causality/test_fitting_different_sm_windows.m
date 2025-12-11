% Double DDM test

clearvars

% Load the table first. We will take advantage of an already existing
% dataset.
threshold_imm = 2; threshold_mob = 2; threshold_pc = 4; id_code = sprintf('imm%d_mob%d_pc%d', threshold_imm, threshold_mob, threshold_pc);
paths = path_generator('folder', 'causality/fitting_windows', 'bouts_id', id_code, 'imfirst', false);
bouts = importdata(fullfile(paths.dataset, 'bouts.mat'));
motion_cache = importdata(fullfile(paths.cache_path, 'motion_cache.mat'));

bouts_proc = data_parser_new(bouts, 'type', 'immobility', 'period', 'loom', 'window', 'le', 'nloom', 2:20);
points.censoring = 10.5;
points.truncation = 0.3;

%  Now we added our vector column to the bouts table

edge_l = -630;
edge_r = 240;
ts = edge_l:(edge_r - 1);
total_length = length(ts);

windows.anchor = 'freeze_onset';
windows.reference = 'fixed_length';
windows.points = [-600 -300 -60 -30 0 5 10 15 20 25 30 35 180];
windows.length = 30;

chunk_len = points.censoring * 60 - 1;
sm_window = nan(height(bouts_proc), total_length);
sm_during_ILI = cell(1, height(bouts_proc));

for idx_trials = 1:height(bouts_proc)

    if strcmp(windows.anchor, 'loom_onset')
        ons = bouts_proc.loom_ts(idx_trials);
    
    elseif strcmp(windows.anchor, 'freeze_onset')
        ons = bouts_proc.onsets(idx_trials);
    end

    off = bouts_proc.ends(idx_trials);
    sum_motion = motion_cache(bouts_proc.fly(idx_trials));

    sm_during_ILI{idx_trials} = sum_motion(ons : ons + chunk_len) ./ 10;
    sm_window(idx_trials, :) = sum_motion(ons + ts) ./ 10;

end

model = 'dddm2';
paths = path_generator('folder', sprintf('causality/fitting_windows/%s', model), 'bouts_id', id_code, 'imfirst', false);
create_output_dirs(paths);

for idx_backframes = windows.points

    backframes = abs(edge_l - idx_backframes) + 1;

    mkdir(fullfile(paths.fig, num2str(backframes)));
    mkdir(fullfile(paths.results, num2str(backframes)));
    paths_backframe.results = fullfile(paths.results, num2str(backframes));
    paths_backframe.fig = fullfile(paths.fig, num2str(backframes));

    if strcmp(windows.reference, 'fixed_length')
        bouts_proc.avg_sm_pre_norm = mean(sm_window(:, backframes:backframes + windows.length - 1), 2);

    elseif strcmp(windows.reference, 'relative')

    end

    %bouts_proc.avg_sm_pre_norm = mean(sm_ts(:, end - backframes:end), 2);
    soc_mot_array = cell2mat(sm_during_ILI)';
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