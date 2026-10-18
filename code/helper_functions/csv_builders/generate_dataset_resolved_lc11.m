% Builds the LC11 dataset_resolved straight from the raw per-fly files of
% 006--lc11/dataset_using_raw_files, the same way generate_dataset_resolved does
% it for CS: one csv per recording, with the surrounding-fly columns (speed,
% angle, motion, their sums and the distances) computed from the raw positions
% and the focal labels recomputed. Each recording is rebuilt in resolve_experiment.
% Raw fly1 is the focal fly, fly2-fly5 are the surrounding ones. The raw files
% have no jumps, so that column is NaN everywhere (no jump, no jump_start).
%
% Files are numbered in the order load_flies_new visits them, by the number of
% moving (NorpA) flies and then by name, and named <recording>.fly<N>.csv. N is
% the fly ID of the bouts table and of the caches. All of it is 50 cm.

clearvars

%% Paths and parameters

raw_dir = '/Users/marcocolnaghi/experimental_data/006--lc11/dataset_using_raw_files';
old_dir = '/Users/marcocolnaghi/experimental_data/006--lc11/dataset_3';
out_dir = '/Users/marcocolnaghi/experimental_data/006--lc11/dataset_resolved';

if ~exist(out_dir, 'dir')
    mkdir(out_dir);
end

% Same as generate_dataset_resolved, see there
params.speed_cap = 75;
params.fly_radius = 1.5;
params.dist_scale = 25;
params.still_vel = 20;
params.freeze_gap = 3;
params.bout_min = 30;

% The raw loom channel has stray pulses (23 of the 1000 recordings), see
% clean_loom_channel in resolve_experiment. The window is anchored on the first
% loom, row 18000, and a loom is 30 frames, 23-31 in the files
params.n_looms = 20;
params.loom_first_row = 18000;
params.loom_bridge = 3;
params.loom_min = 20;
params.loom_max = 40;

%% Experiment list, numbered as load_flies_new visits them

n_moving = 0:4;
token = '*_%dTrh%dNorpA-*.fly1_scored_frames.csv';

experiments = {};
for idx_moving = 1:length(n_moving)
    listing = dir(fullfile(raw_dir, sprintf(token, 4 - n_moving(idx_moving), n_moving(idx_moving))));
    experiments = [experiments, erase(sort({listing.name}), '.fly1_scored_frames.csv')]; %#ok<AGROW>
end

n_exp = numel(experiments);
n_raw = numel(dir(fullfile(raw_dir, '*.fly1_scored_frames.csv')));

assert(n_exp == n_raw, 'generate_dataset_resolved_lc11:missedExperiments', ...
    'The listing found %d experiments but the raw folder holds %d.', n_exp, n_raw);

%% Output columns

sur_cols = cell(1, 12);
for idx_sur = 1:4
    sur_cols(3*idx_sur - 2 : 3*idx_sur) = {sprintf('speed_sur_fly_%d', idx_sur), ...
        sprintf('angle_sur_fly_%d', idx_sur), sprintf('motion_sur_fly_%d', idx_sur)};
end

out_cols = [{'frame', 'looming_bout', 'velocity', 'pixelchange', 'walk_bout', 'freeze_bout', ...
    'jumps', 'low_vel_bout', 'fly_x_mm', 'fly_y_mm'}, sur_cols, ...
    {'sum_speed', 'sum_angle', 'sum_motion', 'freeze', 'dist_1', 'dist_2', 'dist_3', 'dist_4', 'dist_min'}];
n_cols = numel(out_cols);

header = strjoin(out_cols, ',');

% 15 digits and NaN, 7 for the columns derived from the positions, as in CS
is_geometric = startsWith(out_cols, {'angle_', 'motion_', 'dist_'}) | ismember(out_cols, {'sum_angle', 'sum_motion'});
col_fmt = repmat({'%.15g'}, 1, n_cols);
col_fmt(is_geometric) = {'%.7g'};
fmt = [strjoin(col_fmt, ','), '\n'];

%% Resolve and write every experiment

tic
parfor idx_exp = 1:n_exp
    M = resolve_experiment(raw_dir, experiments{idx_exp}, params);

    assert(size(M, 2) == n_cols, 'generate_dataset_resolved_lc11:columnCount', ...
        '%s came out with %d columns instead of %d.', experiments{idx_exp}, size(M, 2), n_cols);

    fid = fopen(fullfile(out_dir, sprintf('%s.fly%d.csv', experiments{idx_exp}, idx_exp)), 'w');
    fprintf(fid, '%s\n', header);
    fprintf(fid, fmt, M');
    fclose(fid);

    fprintf('Fly %d out of %d. \n', idx_exp, n_exp)
end
fprintf('Done in %.1f min. \n', toc / 60)

%% Sanity check

fprintf('%d files in dataset_resolved, %d expected. \n', numel(dir(fullfile(out_dir, '*.csv'))), n_exp)

% The recordings that dataset_3 (Mirjam's files) shares with this set, matched on
% genotype, camera and time stamp. The two come from different tracking, so the
% columns agree closely rather than exactly. The surrounding flies are in a
% different order in the two sets, so only their sums are compared
key = @(name) strjoin(regexp(name, '(M[CH]-\w+?TNT)_.*(Cam[AB]-\d{4}-\d\d-\d\dT\d\d_\d\d_\d\d)', 'tokens', 'once'), ' ');

old_names = {dir(fullfile(old_dir, '*.csv')).name};
old_keys = cellfun(key, old_names, 'UniformOutput', false);
new_names = arrayfun(@(idx) sprintf('%s.fly%d.csv', experiments{idx}, idx), 1:n_exp, 'UniformOutput', false);
new_keys = cellfun(key, new_names, 'UniformOutput', false);

[shared_keys, idx_new, idx_old] = intersect(new_keys, old_keys);
fprintf('%d recordings in common with dataset_3, %d only here, %d only there. \n', ...
    numel(shared_keys), n_exp - numel(shared_keys), numel(old_names) - numel(shared_keys))

% The cleaned loom channel against dataset_3, which has exactly 20 looms everywhere,
% for every shared recording (dataset_3 has one row less, the first 36000 are compared)
loom_equal = false(numel(shared_keys), 1);
opts = detectImportOptions(fullfile(out_dir, new_names{1}));
opts.SelectedVariableNames = {'looming_bout'};
for idx_check = 1:numel(shared_keys)
    new_loom = readmatrix(fullfile(out_dir, new_names{idx_new(idx_check)}), opts);
    old_loom = readmatrix(fullfile(old_dir, old_names{idx_old(idx_check)}), opts);
    n_rows = min(numel(new_loom), numel(old_loom));
    loom_equal(idx_check) = isequal(new_loom(1:n_rows), old_loom(1:n_rows));
end
fprintf('Loom channel identical to dataset_3 in %d out of %d shared recordings. \n', sum(loom_equal), numel(loom_equal))
disp(new_names(idx_new(~loom_equal)))

check_cols = {'velocity', 'sum_speed', 'sum_angle', 'sum_motion'};
for idx_check = round(linspace(1, numel(shared_keys), 6))
    new = readtable(fullfile(out_dir, new_names{idx_new(idx_check)}), 'VariableNamingRule', 'preserve');
    old = readtable(fullfile(old_dir, old_names{idx_old(idx_check)}), 'VariableNamingRule', 'preserve');

    % The raw window is -300 s .. +300 s inclusive, 36001 rows, one more than dataset_3
    n_rows = min(height(new), height(old));
    new = new(1:n_rows, :);
    old = old(1:n_rows, :);

    r = cellfun(@(col) corr(new.(col), old.(col), 'Rows', 'complete'), check_cols);
    fprintf('Fly %d: loom frames equal %d, pixelchange<4 agreement %.4f, correlation with dataset_3 %s. \n', ...
        idx_new(idx_check), isequal(new.looming_bout, old.looming_bout), ...
        mean((new.pixelchange < 4) == (old.pixelchange < 4)), ...
        strjoin(cellfun(@(col, v) sprintf('%s %.3f', col, v), check_cols, num2cell(r), 'UniformOutput', false), ', '))
end
