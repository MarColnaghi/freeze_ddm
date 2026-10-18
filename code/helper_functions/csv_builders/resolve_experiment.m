function M = resolve_experiment(raw_dir, experiment, params)
% Rebuilds one experiment of dataset_resolved from its five raw files. Raw
% fly1 is the focal fly, raw fly2-fly5 are surrounding flies 1-4. Columns
% come out in the order of out_cols in generate_dataset_resolved.

raw = cell(1, 5);
for idx_fly = 1:5
    raw{idx_fly} = read_raw(fullfile(raw_dir, sprintf('%s.fly%d_scored_frames.csv', experiment, idx_fly)));
end

focal = raw{1};
n_frames = numel(focal.velocity);

assert(all(cellfun(@(fly) numel(fly.velocity), raw) == n_frames), 'resolve_experiment:frameCount', ...
    'The raw files of %s do not have the same number of frames.', experiment);

%% Surrounding flies

speed_sur = nan(n_frames, 4);
angle_sur = nan(n_frames, 4);
sur_block = nan(n_frames, 12);

for idx_sur = 1:4
    sur = raw{idx_sur + 1};
    dist_mm = hypot(sur.fly_x_mm - focal.fly_x_mm, sur.fly_y_mm - focal.fly_y_mm);

    % 'includenan' keeps the NaN velocity of frame 0, plain min would turn it into the cap
    speed_sur(:, idx_sur) = min(sur.velocity, params.speed_cap, 'includenan');
    angle_sur(:, idx_sur) = 2 * atan(params.fly_radius ./ dist_mm);

    sur_block(:, 3*idx_sur - 2 : 3*idx_sur) = [speed_sur(:, idx_sur), angle_sur(:, idx_sur), ...
        speed_sur(:, idx_sur) .* angle_sur(:, idx_sur)];
end

% Same formula, and so the same units, as build_distance_cache: mm x 50/3
dist = params.dist_scale ./ tan(angle_sur / 2);

%% Focal fly labels

% A freeze ends at the first change in pixelchange or in velocity. pixelchange
% alone misses most jumps, it often drops to 0 while the fly is in the air.
% Fast frames are removed after the gap filling, so a gap fill never bridges them
fast = focal.velocity >= params.still_vel;

freeze = ~bwareaopen(~(focal.pixelchange == 0), params.freeze_gap + 1) & ~fast; % fill short moving gaps

% freeze_bout is the freeze runs of at least bout_min frames, which is the
% immobility definition of the bouts table at fill-in freeze_gap + 1 (imm3 for
% freeze_gap = 3), bouts of 30 frames or more. pixelchange == 0 stands in for
% the table's pixelchange < thresholds.pc, the same thing while pixelchange is
% only ever 0 or >= 7. If the table's fill-in changes, so must freeze_gap
freeze_bout = bwareaopen(freeze, params.bout_min);

% Upstream the two labels never overlap. Without this the flies repaired by
% the new freeze_bout would be freezing and low-velocity at the same time
low_vel_bout = focal.low_vel_bout;
low_vel_bout(freeze_bout) = 0;

% Only the LC11 raw files need it, which set loom_first_row
looming_bout = focal.looming_bout;
if isfield(params, 'loom_first_row')
    looming_bout = clean_loom_channel(looming_bout, params, experiment);
end

M = [focal.row_index, looming_bout, focal.velocity, focal.pixelchange, focal.walk_bout, ...
    freeze_bout, focal.jumps, low_vel_bout, focal.fly_x_mm, focal.fly_y_mm, sur_block, ...
    sum(speed_sur, 2, 'omitnan'), sum(angle_sur, 2, 'omitnan'), sum(sur_block(:, 3:3:12), 2, 'omitnan'), ...
    freeze, dist, min(dist, [], 2)];

end


function loom = clean_loom_channel(loom, params, experiment)
% The loom channel of the LC11 raw files carries stray pulses: blips of 1-10
% frames, one before the first loom, two long ones of 57 and 105 frames, and in
% one recording a dropped frame splits a loom in two. Gaps of up to loom_bridge
% frames are closed first. Then only the runs of a loom are kept: starting at or
% after loom_first_row (a row index of the file, so 1 more in MATLAB), the first
% loom of the window, and lasting loom_min to loom_max frames (real looms are
% 23-31). What is left must be the n_looms looms.

d = diff([0; double(loom(:) > 0); 0]);
starts = find(d == 1);
ends = find(d == -1) - 1;

join = (starts(2:end) - ends(1:end-1) - 1) <= params.loom_bridge;
starts([false; join]) = [];
ends(join) = [];

len = ends - starts + 1;
is_loom = starts > params.loom_first_row & len >= params.loom_min & len <= params.loom_max;

assert(nnz(is_loom) == params.n_looms, 'resolve_experiment:loomCount', ...
    '%s has %d looms after cleaning the loom channel, expected %d.', experiment, nnz(is_loom), params.n_looms);

loom = zeros(size(loom));
for idx_loom = find(is_loom)'
    loom(starts(idx_loom) : ends(idx_loom)) = 1;
end

end


function fly = read_raw(file_path)
% One field per column, looked up by name since a few raw files carry an
% extra looming_flip column. textscan is much faster than readtable here

fid = fopen(file_path);
assert(fid > 0, 'resolve_experiment:missingFile', 'Cannot open %s.', file_path);

% strtrim drops the \r of the Windows line endings
col_names = strsplit(strtrim(fgetl(fid)), ',');
data = textscan(fid, repmat('%f', 1, numel(col_names)), 'Delimiter', ',', 'EmptyValue', NaN, 'CollectOutput', true);
fclose(fid);

for idx_col = 1:numel(col_names)
    fly.(matlab.lang.makeValidName(col_names{idx_col})) = data{1}(:, idx_col);
end

% The first column is the row index, 'Unnamed: 0' in the CS files and unnamed in the LC11 ones
fly.row_index = data{1}(:, 1);

end
