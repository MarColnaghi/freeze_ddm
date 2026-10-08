% Builds dataset_resolved straight from the raw per-fly files: one csv per
% experiment with the dataset_0 schema, minus the freeze_bout_N columns, plus
% the coordinates of the focal fly. The first column is renamed to frame.
% freeze_bout is recomputed for every fly, which fixes the upstream bug that
% left it empty whenever the fly was still on both the first and last frame.
% freeze and freeze_bout end at the first change in pixelchange or velocity.
% Distances all come from the raw positions, in the dataset_0 units, and the
% columns derived from the positions are written with 7 significant digits.
% Each experiment is rebuilt in resolve_experiment.

clearvars

%% Paths and parameters

raw_dir = '/Users/marcocolnaghi/experimental_data/004--social_ddm/dataset_using_raw_files';
d0_dir  = '/Users/marcocolnaghi/experimental_data/004--social_ddm/dataset_0';
out_dir = '/Users/marcocolnaghi/experimental_data/004--social_ddm/dataset_resolved';

if ~exist(out_dir, 'dir')
    mkdir(out_dir);
end

params.speed_cap = 75;      % surrounding fly speed is clipped here
params.fly_radius = 1.5;    % mm, each fly counts as a 3 mm object in angle_sur_fly
params.dist_scale = 25;     % dist = 25 / tan(angle / 2), the dataset_0 units (mm x 50/3)
params.still_vel = 20;      % mm/s, a frame this fast is moving even with pixelchange == 0
params.freeze_gap = 3;      % freeze fills pixelchange gaps of up to 3 frames
params.bout_gap = 2;        % freeze_bout fills interior pixelchange gaps of up to 2 frames
params.bout_min = 30;       % freeze_bout keeps runs of at least 30 frames

%% Experiment list, numbered as in dataset_0

speeds = [25, 50];
n_moving = 0:4;
token = '1CS%dNorpA%dLC6ChR_5F-%dcm*.fly1_scored_frames.csv';

experiments = {};
for idx_sloom = 1:length(speeds)
    for idx_moving = 1:length(n_moving)
        listing = dir(fullfile(raw_dir, sprintf(token, n_moving(idx_moving), 4 - n_moving(idx_moving), speeds(idx_sloom))));
        experiments = [experiments, erase(sort({listing.name}), '.fly1_scored_frames.csv')]; %#ok<AGROW>
    end
end

n_exp = numel(experiments);
n_raw = numel(dir(fullfile(raw_dir, '*.fly1_scored_frames.csv')));

assert(n_exp == n_raw, 'generate_dataset_resolved:missedExperiments', ...
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

% 15 digits and NaN, as writetable did for dataset_0, except for the columns
% derived from the positions. Those are rounded to 0.01 mm in the raw files,
% so 7 digits are already well below their precision and save ~40% of space
is_geometric = startsWith(out_cols, {'angle_', 'motion_', 'dist_'}) | ismember(out_cols, {'sum_angle', 'sum_motion'});
col_fmt = repmat({'%.15g'}, 1, n_cols);
col_fmt(is_geometric) = {'%.7g'};
fmt = [strjoin(col_fmt, ','), '\n'];

%% Resolve and write every experiment

tic
parfor idx_exp = 1:n_exp
    M = resolve_experiment(raw_dir, experiments{idx_exp}, params);

    assert(size(M, 2) == n_cols, 'generate_dataset_resolved:columnCount', ...
        '%s came out with %d columns instead of %d.', experiments{idx_exp}, size(M, 2), n_cols);

    fid = fopen(fullfile(out_dir, sprintf('%s.fly%d.csv', experiments{idx_exp}, idx_exp)), 'w');
    fprintf(fid, '%s\n', header);
    fprintf(fid, fmt, M');
    fclose(fid);

    fprintf('Fly %d out of %d. \n', idx_exp, n_exp)
end
fprintf('Done in %.1f min. \n', toc / 60)

%% Sanity check against dataset_0

names = arrayfun(@(idx) sprintf('%s.fly%d.csv', experiments{idx}, idx), 1:n_exp, 'UniformOutput', false);

fprintf('%d files in dataset_resolved, %d expected. \n', numel(dir(fullfile(out_dir, '*.csv'))), n_exp)
fprintf('%d out of %d names match dataset_0. \n', sum(cellfun(@(name) isfile(fullfile(d0_dir, name)), names)), n_exp)

% Upstream freeze_bout was empty for these flies, still on the first and last frame
bug_flies = [3, 88, 135, 290, 516, 529, 551, 560, 599, 606, 615, 616, 626, 640, 653, 656, ...
    681, 692, 736, 785, 845, 904, 958];

opts = detectImportOptions(fullfile(out_dir, names{1}));
opts.SelectedVariableNames = {'freeze_bout'};
bug_frames = arrayfun(@(idx) sum(readmatrix(fullfile(out_dir, names{idx}), opts)), bug_flies);
fprintf('freeze_bout is no longer empty in %d out of %d previously bugged flies. \n', sum(bug_frames > 0), numel(bug_flies))

% A few flies in full: every shared column should match dataset_0 to within the
% 7 digit rounding (relative diff <= 5e-7), except the recomputed labels, which
% only differ on the fast frames (and for freeze_bout in the bugged flies), and
% the distances. dataset_0 took those from the aligned pixel coordinates for 370
% flies, which agree with the raw positions to < 0.1 mm, except flies 149, 451
% and 944 where the aligned coordinates are off by several mm
label_cols = {'freeze', 'freeze_bout', 'low_vel_bout'};
dist_cols = {'dist_1', 'dist_2', 'dist_3', 'dist_4', 'dist_min'};

for idx_fly = [bug_flies(1), 149, randperm(n_exp, 5)]
    new = readtable(fullfile(out_dir, names{idx_fly}), 'VariableNamingRule', 'preserve');
    old = readtable(fullfile(d0_dir, names{idx_fly}), 'VariableNamingRule', 'preserve');
    old.Properties.VariableNames{1} = 'frame';

    shared = setdiff(intersect(new.Properties.VariableNames, old.Properties.VariableNames), [label_cols, dist_cols]);
    max_diff = max(cellfun(@(col) max(abs(new.(col) - old.(col)) ./ max(abs(old.(col)), 1e-12), [], 'omitnan'), shared));
    nan_mismatch = sum(cellfun(@(col) sum(isnan(new.(col)) ~= isnan(old.(col))), shared));
    dist_diff = max(cellfun(@(col) max(abs(new.(col) - old.(col))), dist_cols)) / (50/3); % in mm

    fprintf(['Fly %d: max relative diff %.2g over %d shared columns, %d NaN mismatches, distances within %.2g mm, ', ...
        'freeze agreement %.4f, freeze_bout agreement %.4f (%d vs %d frames). \n'], idx_fly, max_diff, numel(shared), ...
        nan_mismatch, dist_diff, mean(new.freeze == old.freeze), mean(new.freeze_bout == old.freeze_bout), ...
        sum(new.freeze_bout), sum(old.freeze_bout))
end
