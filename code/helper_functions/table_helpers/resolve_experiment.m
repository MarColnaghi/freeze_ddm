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
freeze_bout = freeze_bout_mask(focal.pixelchange, fast, params);

% Upstream the two labels never overlap. Without this the flies repaired by
% the new freeze_bout would be freezing and low-velocity at the same time
low_vel_bout = focal.low_vel_bout;
low_vel_bout(freeze_bout) = 0;

M = [focal.Unnamed_0, focal.looming_bout, focal.velocity, focal.pixelchange, focal.walk_bout, ...
    freeze_bout, focal.jumps, low_vel_bout, focal.fly_x_mm, focal.fly_y_mm, sur_block, ...
    sum(speed_sur, 2, 'omitnan'), sum(angle_sur, 2, 'omitnan'), sum(sur_block(:, 3:3:12), 2, 'omitnan'), ...
    freeze, dist, min(dist, [], 2)];

end


function freeze_bout = freeze_bout_mask(pixelchange, fast, params)
% Upstream freeze_bout without its start/end pairing bug, which emptied it
% whenever the fly was still on both the first and the last frame. Frames
% with pixelchange == 0 and interior gaps of up to bout_gap frames filled,
% minus the fast ones, keeping only runs of at least bout_min frames

still = pixelchange == 0;

freeze_bout = false(size(still));
if ~any(still)
    return
end

% Gaps are filled only between two still stretches, never before the first
% still frame or after the last one
gaps = ~bwareaopen(~still, params.bout_gap + 1) & ~still;
gaps([1:find(still, 1) - 1, find(still, 1, 'last') + 1:end]) = false;

freeze_bout = bwareaopen((still | gaps) & ~fast, params.bout_min);

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

end
