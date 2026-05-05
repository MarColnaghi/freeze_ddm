% Define the directories
folder_dataset1 = '/Users/marcocolnaghi/experimental_data/004--social_ddm/dataset_1';
folder_dataset0 = '/Users/marcocolnaghi/experimental_data/004--social_ddm/dataset_0';
folder_aligned  = '/Volumes/Elements/SocialDDM-CS/aligned';

% Create dataset_0 if it doesn't already exist
if ~exist(folder_dataset0, 'dir')
    mkdir(folder_dataset0);
    fprintf('Created new directory: %s\n', folder_dataset0);
end

% Get all csv files from dataset_1
files = dir(fullfile(folder_dataset1, "*.csv"));

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
    current_fly = fly_num(i);
    fprintf('Processing fly %d: %s\n', current_fly, fname);
    
    % Define source (dataset_1) and destination (dataset_0) paths
    src_fpath  = fullfile(folder_dataset1, fname);
    dest_fpath = fullfile(folder_dataset0, fname);
    
    % Read the original table from dataset_1
    T = readtable(src_fpath, 'TextType','string', 'PreserveVariableNames', true);
    
    % Flag to track whether we successfully computed distances
    dist_computed = false;
    
    % ==========================================
    % METHOD 1: Try Coordinate Method First
    % ==========================================
    fly_coord_folder = fullfile(folder_aligned, sprintf('fly%d', current_fly));
    coord_files = dir(fullfile(fly_coord_folder, 'aligned_*_coordinates.csv'));
    
    if ~isempty(coord_files)
        coord_fpath = fullfile(fly_coord_folder, coord_files(1).name);
        T_coords = readtable(coord_fpath);
        
        % Assuming first 10 columns are x1,y1 to x5,y5
        coords_mat = table2array(T_coords(:, 1:10)); 
        
        % Check for row mismatch
        if height(T) == size(coords_mat, 1)
            x1 = coords_mat(:, 1);
            y1 = coords_mat(:, 2);
            
            dist_1 = sqrt((coords_mat(:, 3) - x1).^2 + (coords_mat(:, 4) - y1).^2);
            dist_2 = sqrt((coords_mat(:, 5) - x1).^2 + (coords_mat(:, 6) - y1).^2);
            dist_3 = sqrt((coords_mat(:, 7) - x1).^2 + (coords_mat(:, 8) - y1).^2);
            dist_4 = sqrt((coords_mat(:, 9) - x1).^2 + (coords_mat(:, 10) - y1).^2);
            
            dist_computed = true;
            fprintf('  -> Used Coordinates to compute distance.\n');
        else
            fprintf('  -> Row mismatch (Dataset: %d, Coords: %d). Falling back to Angle method.\n', height(T), size(coords_mat, 1));
        end
    else
        fprintf('  -> No coordinate file found. Falling back to Angle method.\n');
    end
    
    % ==========================================
    % METHOD 2: Fallback to Angle Method
    % ==========================================
    if ~dist_computed
        % Verify the columns exist before trying to use them
        ang_cols = {'angle_sur_fly_1', 'angle_sur_fly_2', 'angle_sur_fly_3', 'angle_sur_fly_4'}; % CHANGE IF NEEDED
        
        if all(ismember(ang_cols, T.Properties.VariableNames))
            dist_1 = 25 ./ tan(T{:, 'angle_sur_fly_1'} / 2);
            dist_2 = 25 ./ tan(T{:, 'angle_sur_fly_2'} / 2);
            dist_3 = 25 ./ tan(T{:, 'angle_sur_fly_3'} / 2);
            dist_4 = 25 ./ tan(T{:, 'angle_sur_fly_4'} / 2);
            
            dist_computed = true;
            fprintf('  -> Used Angular Size to compute distance.\n');
        else
            fprintf('  -> ERROR: Neither coordinates nor angular size columns found! Skipping calculation.\n');
        end
    end
    
    % ==========================================
    % APPEND AND SAVE
    % ==========================================
    if dist_computed
        T.dist_1 = dist_1;
        T.dist_2 = dist_2;
        T.dist_3 = dist_3;
        T.dist_4 = dist_4;
        T.dist_min = min([dist_1, dist_2, dist_3, dist_4], [], 2);
    end
    
    % Save to dataset_0 (saves updated table, or just copies it if both methods failed)
    writetable(T, dest_fpath);
end

fprintf('\nDone processing all files.\n');