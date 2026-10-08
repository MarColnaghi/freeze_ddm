function arenas = build_arena_cache()
% BUILD_ARENA_CACHE - Fits the arena circle of every experiment
%   The arena sits at a different place in the camera frame from one recording
%   to the next. Each experiment's arena is the circle (Kasa fit) through the
%   convex hull of every position of all five flies, focal and surrounding,
%   over the whole recording, read from the raw per-fly files.
%
%   A few tracks are broken, e.g. surrounding flies 2-4 of fly 639 run along
%   axis-aligned lines well outside the arena, and a single one drags the
%   hull. Each fly's own hull is fitted first, and only the flies whose hull
%   is round (spread around their own circle <= max_fly_sd) enter the fit.
%
% Returns:
%   arenas: table, one row per experiment, keyed by the focal fly ID used in
%   bouts. centre_x, centre_y and radius are in mm in the raw tracking frame.
%   hull_sd is the spread of the hull vertices around the fitted circle.
%   fly_hull_sd holds the same spread for each fly's own fit (raw fly1-fly5,
%   fly1 the focal fly) and fly_used which of them entered the arena fit

raw_dir = '/Users/marcocolnaghi/experimental_data/004--social_ddm/dataset_using_raw_files';
resolved_dir = '/Users/marcocolnaghi/experimental_data/004--social_ddm/dataset_resolved';

% dataset_resolved names each experiment <experiment>.fly<N>.csv, with N the
% fly ID that load_flies_new puts in bouts.fly
listing = dir(fullfile(resolved_dir, '*.csv'));
tokens = regexp({listing.name}, '^(.*)\.fly(\d+)\.csv$', 'tokens', 'once');
experiment = cellfun(@(t) t{1}, tokens, 'UniformOutput', false)';
fly = cellfun(@(t) str2double(t{2}), tokens)';

max_fly_sd = 1; % mm, clean tracks sit at ~0.3

n_exp = numel(fly);
centre = nan(n_exp, 2);
radius = nan(n_exp, 1);
hull_sd = nan(n_exp, 1);
fly_hull_sd = nan(n_exp, 5);
fly_used = false(n_exp, 5);

parfor idx_exp = 1:n_exp
    xy = cell(5, 1);
    sd_exp = nan(1, 5);
    for idx_fly = 1:5
        p = read_xy(fullfile(raw_dir, sprintf('%s.fly%d_scored_frames.csv', experiment{idx_exp}, idx_fly)));
        xy{idx_fly} = p(all(~isnan(p), 2), :);
        [~, ~, sd_exp(idx_fly)] = fit_hull_circle(xy{idx_fly});
    end

    used = sd_exp <= max_fly_sd;
    if ~any(used)
        used = sd_exp == min(sd_exp);
    end

    [c, r, sd] = fit_hull_circle(vertcat(xy{used}));
    centre(idx_exp, :) = c;
    radius(idx_exp) = r;
    hull_sd(idx_exp) = sd;
    fly_hull_sd(idx_exp, :) = sd_exp;
    fly_used(idx_exp, :) = used;
end

arenas = table(fly, experiment, centre(:, 1), centre(:, 2), radius, hull_sd, fly_hull_sd, fly_used, ...
    'VariableNames', {'fly', 'experiment', 'centre_x', 'centre_y', 'radius', 'hull_sd', 'fly_hull_sd', 'fly_used'});
arenas = sortrows(arenas, 'fly');

paths = path_generator();
save(fullfile(paths.cache_path, 'arena_cache.mat'), 'arenas')

end


function [c, r, sd] = fit_hull_circle(xy)
% Kasa circle through the convex hull vertices, and their spread around it

h = convhull(xy(:, 1), xy(:, 2));
xh = xy(h, 1);
yh = xy(h, 2);
s = [xh, yh, ones(numel(h), 1)] \ -(xh.^2 + yh.^2);
c = -s(1:2)' / 2;
r = sqrt(sum(c.^2) - s(3));
sd = std(hypot(xh - c(1), yh - c(2)) - r);

end


function xy = read_xy(file_path)
% Only the two position columns, looked up by name since a few raw files carry
% an extra looming_flip column

fid = fopen(file_path);
assert(fid > 0, 'build_arena_cache:missingFile', 'Cannot open %s.', file_path);

% strtrim drops the \r of the Windows line endings
col_names = strsplit(strtrim(fgetl(fid)), ',');
idx_x = find(strcmp(col_names, 'fly_x_mm'));
idx_y = find(strcmp(col_names, 'fly_y_mm'));

fmt = repmat({'%*f'}, 1, numel(col_names));
fmt([idx_x, idx_y]) = {'%f'};
data = textscan(fid, [fmt{:}], 'Delimiter', ',', 'EmptyValue', NaN, 'CollectOutput', true);
fclose(fid);

% textscan returns the columns in file order
xy = data{1};
if idx_y < idx_x
    xy = xy(:, [2 1]);
end

end
