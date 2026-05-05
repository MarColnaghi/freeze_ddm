function caches = build_pc_cache_from_bouts(bouts, varargin)

% Setup Input Parser
opt = inputParser;
addParameter(opt, 'save_external_disk', false);
parse(opt, varargin{:});

% Data Processing Setup
thresholds = define_thresholds;
thresholds.le_window_fl = [5 55];
thresholds.le_window_sl = [5 55];

% Process bouts with thresholds
bouts = bouts_formatting(bouts, thresholds);

% Initialize the Map for all flies
caches = containers.Map('KeyType', 'double', 'ValueType', 'any');

% Define Base Paths
external_base_path = '/Volumes/Elements/SocialDDM-CS/aligned';
local_save_path = '/Users/marcocolnaghi/PhD/freeze_ddm/datasets/caches/freeze_cache_from_bout.mat';

% Main Loop
for idx_fly = unique(bouts.fly)'
    % Filter data for current fly
    bouts_fly = bouts(bouts.fly == idx_fly, :);
    bouts_fly = sortrows(bouts_fly, {'onsets'});
    
    % Calculate durations and handle invalid values
    dur = bouts_fly.durations;
    dur(dur <= 0 | isnan(dur)) = 0;
    
    % Logical processing for 'x' vector
    base_val = double(bouts_fly.type);
    
    % Condition: Long lasting bouts (>= 30) with specific LE and frozen states
    condition = (bouts_fly.le == true & ...
                 bouts_fly.type == true & ...
                 bouts_fly.frozen_start == false & ...
                 bouts_fly.durations >= 30);
    base_val(condition) = 2;
    
    % Create the expanded time-series vector
    x = repelem(base_val, dur);
    x = [x; x(end)]; % Maintain indexing consistency
    
    % Store in the master Map
    caches(idx_fly) = x(:)';
    
    % External Disk Save Logic
    if opt.Results.save_external_disk
        % Construct path to specific fly folder: /Volumes/.../aligned/fly1/
        fly_folder = fullfile(external_base_path, sprintf('fly%d', idx_fly));
        
        % Check if the folder exists on the external drive
        if exist(fly_folder, 'dir')
            % Define filename for this specific fly
            external_file_path = fullfile(fly_folder, 'freeze_cache_from_bout.csv');
            
            % Save ONLY this fly's vector (x) to the external drive
            writematrix(x, external_file_path);
            fprintf('Saved: Fly %d data to external disk.\n', idx_fly);
        else
            fprintf('Warning: Directory not found, skipping external save: %s\n', fly_folder);
        end
    end
end

% Master Local Save
% Check if the local directory exists before saving the full Map
[local_dir, ~, ~] = fileparts(local_save_path);
if ~exist(local_dir, 'dir')
    mkdir(local_dir);
end

save(local_save_path, 'caches');
fprintf('Success: Master cache saved to %s\n', local_save_path);