% Double DDM test - Refined Sliding Window Analysis
clearvars

% 1. Initialization and Data Loading
threshold_imm = 2; threshold_mob = 2; threshold_pc = 4; 
id_code = sprintf('imm%d_mob%d_pc%d', threshold_imm, threshold_mob, threshold_pc);
paths = path_generator('folder', 'causality/fitting_windows', 'bouts_id', id_code, 'imfirst', false);

bouts = importdata(fullfile(paths.dataset, 'bouts.mat'));
motion_cache = importdata(fullfile(paths.cache_path, 'motion_cache.mat'));

% 2. Behavioral Formatting
thresholds = define_thresholds('le_window', struct('le_window_sl', [5 55], 'le_window_fl', [5 55]));
bouts = bouts_formatting(bouts, thresholds);
bouts_proc = data_parser_new(bouts, 'type', 'immobility', 'period', 'loom', 'window', 'le', 'nloom', 1:20, 'min_dur', 30);

% 3. Window and Timeline Configuration
points.censoring = 10.5;
points.truncation = min(bouts_proc.durations_s);

edge_l = -630;
edge_r = 630; % Increased to accommodate windows.points up to 600
ts = edge_l:(edge_r - 1);

% Your points will now all be within bounds
windows.points = [-600 -300 -60 -30 -10 10 20 30 40 50 60 120 180 240 360 600];
onset_idx = find(ts == 0); % The reference "Zero" point

windows.anchor = 'freeze_onset';
windows.reference = 'fixed_length'; % or 'cumulative'
windows.length = '30'; % Set a numeric string for fixed_length (e.g., '30' frames)

% 4. Pre-calculating Motion Matrix (sm_window)
total_length = length(ts);
chunk_len = points.censoring * 60 - 1;
sm_window = nan(height(bouts_proc), total_length);
sm_during_ILI = cell(1, height(bouts_proc));

for idx_trials = 1:height(bouts_proc)
    if strcmp(windows.anchor, 'loom_onset')
        ons = bouts_proc.loom_ts(idx_trials);
    elseif strcmp(windows.anchor, 'freeze_onset')
        ons = bouts_proc.onsets(idx_trials);
    end
    
    sum_motion = motion_cache(bouts_proc.fly(idx_trials));
    
    % Store the chunk for the DDM fitting
    sm_during_ILI{idx_trials} = sum_motion(ons : ons + chunk_len) ./ 10;
    
    % Store the full window for the sliding regressor
    sm_window(idx_trials, :) = sum_motion(ons + ts) ./ 10;
end

% 5. Directory Setup
model = 'dddm2';
paths = path_generator('folder', sprintf('causality/fitting_windows/%s/%s', model, windows.reference), 'bouts_id', id_code, 'imfirst', false);
create_output_dirs(paths, windows);

% 6. Main Iteration Loop
for idx_point = windows.points
    % Find the index in 'ts' corresponding to our current time point
    target_idx = find(ts == idx_point);
    
    if isempty(target_idx)
        fprintf('\n Skipping point %d: Out of bounds.', idx_point);
        continue;
    end

    % --- Calculate the Regressor (smp) ---
    if strcmp(windows.reference, 'fixed_length')
        f_len = str2double(windows.length);
        % Ensure we don't exceed matrix bounds
        end_idx = min(target_idx + f_len - 1, size(sm_window, 2));
        bouts_proc.smp = mean(sm_window(:, target_idx : end_idx), 2);
        fprintf('\n Mode: FIXED | Window: ts [%d to %d]', ts(target_idx), ts(end_idx));

    elseif strcmp(windows.reference, 'cumulative')
        if idx_point < 0
            % From the past up to the onset
            bouts_proc.smp = mean(sm_window(:, target_idx : onset_idx), 2);
            fprintf('\n Mode: CUMULATIVE (PRE) | Window: ts [%d to 0]', idx_point);
        else
            % From onset up to the future point
            bouts_proc.smp = mean(sm_window(:, onset_idx : target_idx), 2);
            fprintf('\n Mode: CUMULATIVE (POST) | Window: ts [0 to %d]', idx_point);
        end
    end

    % 7. Run Fitting
    paths_backframe.results = fullfile(paths.results, num2str(idx_point));
    paths_backframe.fig = fullfile(paths.fig, num2str(idx_point));
    if ~exist(paths_backframe.results, 'dir'), mkdir(paths_backframe.results); end
    if ~exist(paths_backframe.fig, 'dir'), mkdir(paths_backframe.fig); end

    soc_mot_array = cell2mat(sm_during_ILI)';
    extra.soc_mot_array = soc_mot_array;
    
    model_results = run_fitting_newer(bouts_proc, points, model, paths, 'export', false,...
        'bads_display', 'none', 'pass_ndt', false, 'n_bads', 3, 'extra', extra);
    
    % 8. Save Results
    model_results.bouts_path = paths_backframe.results; 
    model_results.paths = paths.fig;
    model_results.motion_cache_path = fullfile(paths.cache_path, 'motion_cache.mat');
    
    save(fullfile(paths_backframe.results, sprintf('fit_results_%s.mat', model)), '-struct', 'model_results');
    save(fullfile(paths_backframe.results, 'extra.mat'), 'extra');
    save(fullfile(paths_backframe.results, 'surrogate.mat'), 'bouts_proc');
end

function create_output_dirs(paths, windows)
    if ~exist(paths.fig, 'dir'), mkdir(paths.fig); end
    if ~exist(paths.results, 'dir'), mkdir(paths.results); end
    
    run_folders = dir(fullfile(paths.results, 'run*'));
    
    if isempty(run_folders)
        next_run = 1;
    else
        % Extract digits, convert to double, and find the max
        matches = regexp({run_folders.name}, '\d+', 'match');
        % Filter out empty matches and convert the first match of each folder to a number
        run_nums = cellfun(@(x) str2double(x{1}), matches(~cellfun(@isempty, matches)));
        next_run = max(run_nums) + 1;
    end
    
    run_name = sprintf('run%02d_ref-%s_anchor-%s', next_run, windows.reference, windows.anchor);
    paths.results = fullfile(paths.results, run_name);
    paths.fig = fullfile(paths.fig, run_name);
    
    mkdir(paths.results); mkdir(paths.fig);
    assignin('caller', 'paths', paths);
end

function val = if_else(condition, true_val, false_val)
    if condition, val = true_val; else, val = false_val; end
end