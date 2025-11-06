function caches = build_caches(var_name)
% BUILD_SOCIAL_SIGNAL_MAP - Minimal version
%   sigMap = build_social_signal_map("/Users/marcocolnaghi/experimental_data/004--social_ddm/dataset_1")
%   sigMap = build_social_signal_map(folder, "speed_sur_fly_2")  % choose another column if you want
%
% - folder: folder with CSVs (no recursion)
% - socialCol: column to read (default "speed_sur_fly_1")
%
% Returns:
%   sigMap: containers.Map where key = filename (with extension), value = numeric vector

folder = '/Users/marcocolnaghi/experimental_data/004--social_ddm/dataset_1';

files = dir(fullfile(folder, "*.csv"));
caches = containers.Map('KeyType','double','ValueType','any');

fly_num = nan(numel(files), 1);

for i = 1:numel(files)
    fname = files(i).name;

    fly_num_str = regexp(fname, '(?<=fly)\d+', 'match', 'once');
    fly_num(i) = str2double(fly_num_str);

end

[~, idx] = sort(fly_num);

for i = idx'

    fprintf('fly %d \n', fly_num(i));
    fname = files(i).name;
    fname
    fpath = fullfile(files(i).folder, fname);

    T = readtable(fpath, 'TextType','string', 'PreserveVariableNames', true);

    x = T.(var_name);
    % Coerce to numeric if needed
    if isstring(x) || ischar(x) || iscellstr(x)
        x = str2double(string(x));
    else
        x = double(x);
    end

    caches(fly_num(i)) = x(:);

end

if strcmp(var_name, 'sum_motion')
    save('/Users/marcocolnaghi/PhD/freeze_ddm/datasets/caches/motion_cache.mat', 'caches')
elseif strcmp(var_name, 'pixelchange')
    save('/Users/marcocolnaghi/PhD/freeze_ddm/datasets/caches/pixel_cache.mat', 'caches')
elseif strcmp(var_name, 'looming_bout')
    save('/Users/marcocolnaghi/PhD/freeze_ddm/datasets/caches/loom_cache.mat', 'caches')

end