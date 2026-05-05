function caches = build_caches(var_name)
% BUILD_CACHES - Builds caches for signal data from fly CSVs
%   - folder: folder with CSVs
%   - var_name: column name (or 'sur_speed'/'sur_angle' for multi-column)
%
% Returns:
%   caches: containers.Map where key = fly number, value = numeric vector or matrix

folder = '/Users/marcocolnaghi/experimental_data/004--social_ddm/dataset_0';
files = dir(fullfile(folder, "*.csv"));

% Initialize the map
caches = containers.Map('KeyType','double','ValueType','any');

% Extract fly numbers for sorting
fly_num = nan(numel(files), 1);
for i = 1:numel(files)
    fname = files(i).name;
    fly_num_str = regexp(fname, '(?<=fly)\d+', 'match', 'once');
    fly_num(i) = str2double(fly_num_str);
end

% Sort indices by fly number
[~, idx] = sort(fly_num);

for i = idx'
    fname = files(i).name;
    fprintf('Processing fly %d: %s\n', fly_num(i), fname);

    fpath = fullfile(files(i).folder, fname);
    T = readtable(fpath, 'TextType','string', 'PreserveVariableNames', true);

    % --- Handle Multi-column extraction for 'sur' cases ---
    if strcmp(var_name, 'sur_speed')
        col_names = {'speed_sur_fly_1', 'speed_sur_fly_2', 'speed_sur_fly_3', 'speed_sur_fly_4'};
        x = T{:, col_names};

    elseif strcmp(var_name, 'sur_distance')
        col_names = {'dist_1', 'dist_2', 'dist_3', 'dist_4'};

        % Check if all specified columns exist in table T
        if all(ismember(col_names, T.Properties.VariableNames))
            x = T{:, col_names};
        else
            % If they are missing, create an array of NaNs matching T's height
            x = nan(height(T), numel(col_names));
        end

    elseif strcmp(var_name, 'sur_angle')
        col_names = {'angle_sur_fly_1', 'angle_sur_fly_2', 'angle_sur_fly_3', 'angle_sur_fly_4'};
        x = T{:, col_names};
    else
        % Default single column case
        % Check if all specified columns exist in table T
        if all(ismember(var_name, T.Properties.VariableNames))
            x = T{:, var_name};
        else
            % If they are missing, create an array of NaNs matching T's height
            x = nan(height(T), numel(var_name));
        end
    end

    % Coerce to numeric (handles strings or numeric columns)
    if isstring(x) || ischar(x) || iscell(x)
        x = str2double(string(x));
    else
        x = double(x);
    end

    % Store in map: use (:) for vectors to ensure column orientation,
    % but keep as matrix (N x 4) for surrounding fly data.
    if size(x, 2) > 1
        caches(fly_num(i)) = x;
    else
        caches(fly_num(i)) = x(:);
    end
end

% --- Save Logic ---
base_path = '/Users/marcocolnaghi/PhD/freeze_ddm/datasets/caches/';
switch var_name
    case 'sum_motion'
        save(fullfile(base_path, 'motion_cache.mat'), 'caches')
    case 'sum_speed'
        save(fullfile(base_path, 'sumspeed_cache.mat'), 'caches')
    case 'sum_angle'
        save(fullfile(base_path, 'sumangle_cache.mat'), 'caches')
    case 'dist_min'
        save(fullfile(base_path, 'mindist_cache.mat'), 'caches')
    case 'pixelchange'
        save(fullfile(base_path, 'pixel_cache.mat'), 'caches')
    case 'looming_bout'
        save(fullfile(base_path, 'loom_cache.mat'), 'caches')
    case 'velocity'
        save(fullfile(base_path, 'speed_cache.mat'), 'caches')
    case 'freeze_bout'
        save(fullfile(base_path, 'freeze_cache.mat'), 'caches')
    case 'sur_speed'
        save(fullfile(base_path, 'surspeed_cache.mat'), 'caches')
    case 'sur_angle'
        save(fullfile(base_path, 'surangle_cache.mat'), 'caches')
    case 'sur_distance'
        save(fullfile(base_path, 'surdistance_cache.mat'), 'caches')
end

end