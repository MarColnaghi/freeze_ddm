
clearvars

% Load the table first. We will take advantage of an already existing
% dataset.
threshold_imm = 2; threshold_mob = 2; threshold_pc = 4; id_code = sprintf('imm%d_mob%d_pc%d', threshold_imm, threshold_mob, threshold_pc);
paths = path_generator('folder', 'social_motion/clustering', 'bouts_id', id_code, 'imfirst', false);
bouts = importdata(fullfile(paths.dataset, 'bouts.mat'));

motion_cache = importdata(fullfile(paths.cache_path, 'motion_cache.mat'));
speed_cache = importdata(fullfile(paths.cache_path, 'speed_cache.mat'));
angle_cache = importdata(fullfile(paths.cache_path, 'angle_cache.mat'));
loom_cache = importdata(fullfile(paths.cache_path, 'loom_cache.mat'));

%% 
idx_fly = 425;
loom_ts = find(diff(loom_cache(idx_fly)) == 1) + 1;

file_pattern = sprintf('*.fly%d*.csv', idx_fly); 
full_pattern = fullfile(paths.exp_data, file_pattern);
file_info = dir(full_pattern);

if ~isempty(file_info)
    actual_filename = file_info(1).name; 
    data = readtable(fullfile(paths.exp_data, actual_filename));
else
    warning('No file found for fly %d', idx_fly);
end

%

a = corrcoef(speed_cache(idx_fly), motion_cache(idx_fly))
b = corrcoef(speed_cache(idx_fly), angle_cache(idx_fly))
c = corrcoef(angle_cache(idx_fly), motion_cache(idx_fly))


fh = figure('color','w','Position',[100, 100, 950, 1200]);
t = tiledlayout(3, 1, 'TileSpacing', 'tight');
ax = gobjects(3, 1); % Pre-allocate axes handles

% Define labels and limits for each tile to keep the loop clean
field_types = {'speed', 'angle', 'motion'};
y_limits = {[0 150], [0 3], [0 50]};

dt = 1/60;
t = dt:dt:600;

for i = 1:3
    ax(i) = nexttile;
    hold on
    
    current_type = field_types{i};
    
    if i < 3
        % Logic for Speed (1) and Angle (2)
        for idx_surr_flies = 1:4
            surr_data = data.(sprintf('%s_sur_fly_%d', current_type, idx_surr_flies));
            plot(t, surr_data, 'LineWidth', 0.6);
        end
        
        % Plot the main fly data (from cache)
        main_data = eval(sprintf('%s_cache(idx_fly)', current_type)); 
        plot(t, main_data, 'k-', 'LineWidth', 1.2);
        
    else
        % Specific logic for the Motion tile (3)
        % Calculate surrogate motion
        fly_motion = [];
        for idx_surr_flies = 1:4
            s = data.(sprintf('speed_sur_fly_%d', idx_surr_flies));
            a = data.(sprintf('angle_sur_fly_%d', idx_surr_flies));
            fly_motion(idx_surr_flies, :) = s .* a;
        end
        
        plot(t, sum(fly_motion, 1), 'r-', 'LineWidth', 1.2);
    end
    
    xline(loom_ts ./ 60, 'k-.', 'LineWidth', 1)

    % Apply formatting
    x_max = size(speed_cache(idx_fly), 1); % Shared X limit
    
    apply_generic(gca, 'ylim', y_limits{i}, 'xlim', [300 600], 'font_size', 20, 'tick_length', 0.015, 'line_width', 2);
    ylabel(field_types{i})
end

linkaxes(ax, 'x');

