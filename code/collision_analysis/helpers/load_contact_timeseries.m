function data = load_contact_timeseries(ls)
% LOAD_CONTACT_TIMESERIES  Shared data loading for collision/contact analyses.
%
%   data = load_contact_timeseries(ls) loads the cached fly-by-frame
%   timeseries for every fly recorded at loom speed LS, together with the
%   per-fly loom-onset frames and the number of moving conspecifics. Used by
%   the scripts in collision_analysis/contacts/ so the loading boilerplate
%   lives in one place.
%
%   Fields of DATA:
%     sm, fs, fr, pc, md, loom : fly x frame matrices (motion, speed, freeze,
%                                pixel-change, min-distance, loom indicator)
%     selected_flies           : global fly ids included (column, sorted)
%     n_moving_flies           : # moving conspecifics per selected fly
%     loom_times               : nFlies x nLooms onset frame of each loom
%     n_looms                  : looms per fly (asserted equal across flies)
%     n_flies                  : number of selected flies

% Paths + cached timeseries
paths = path_generator('folder', fullfile('descriptive', 'rasters'));

sm_cache   = importdata(fullfile(paths.cache_path, 'motion_cache.mat'));
pc_cache   = importdata(fullfile(paths.cache_path, 'pixel_cache.mat'));
fr_cache   = importdata(fullfile(paths.cache_path, 'freeze_cache.mat'));
fs_cache   = importdata(fullfile(paths.cache_path, 'speed_cache.mat'));
md_cache   = importdata(fullfile(paths.cache_path, 'mindist_cache.mat'));
loom_cache = importdata(fullfile(paths.cache_path, 'loom_cache.mat'));

% Bouts
thresholds = define_thresholds;
bouts      = importdata(fullfile(paths.dataset, 'bouts.mat'));
bouts      = bouts_formatting(bouts, thresholds);

% Flies recorded at the requested loom speed
selected_flies = unique(bouts.fly(bouts.sloom == ls, :));
n_moving_flies = accumarray(bouts.fly, bouts.moving_flies, [], @(x) x(1));
n_moving_flies = n_moving_flies(selected_flies);

% Fly x frame matrices
sm_mat   = cache2mat(sm_cache,   'selected_flies', selected_flies');
fs_mat   = cache2mat(fs_cache,   'selected_flies', selected_flies');
fr_mat   = cache2mat(fr_cache,   'selected_flies', selected_flies');
pc_mat   = cache2mat(pc_cache,   'selected_flies', selected_flies');
md_mat   = cache2mat(md_cache,   'selected_flies', selected_flies');
loom_mat = cache2mat(loom_cache, 'selected_flies', selected_flies');

% Loom onset frames (one row per fly)
loom_ts         = diff(loom_mat, [], 2) == 1;
n_flies         = size(loom_mat, 1);
n_looms_per_fly = sum(loom_ts, 2);
n_looms         = n_looms_per_fly(1);
assert(all(n_looms_per_fly == n_looms), 'Flies do not all have the same number of looms.');

loom_times = nan(n_flies, n_looms);
for f = 1:n_flies
    loom_times(f, :) = find(loom_ts(f, :));
end

% Package
data = struct( ...
    'sm', sm_mat, 'fs', fs_mat, 'fr', fr_mat, 'pc', pc_mat, 'md', md_mat, 'loom', loom_mat, ...
    'selected_flies', selected_flies, 'n_moving_flies', n_moving_flies, ...
    'loom_times', loom_times, 'n_looms', n_looms, 'n_flies', n_flies);
end
