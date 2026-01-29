%% Random Effects Correlation Matrix Heatmap

disp('Generating random effects correlation heatmap...');

% =========================================================================
% 1. CONFIGURATION (can be ignored)
%    Get the table
% =========================================================================
id_code = 'imm3_mob3_pc4';
paths_out = path_generator('folder', 'descriptive/fd_durs', 'bouts_id', id_code, 'imfirst', false);

% col = cmapper();
thresholds = define_thresholds;
thresholds.le_window_fl = [5 40];
thresholds.le_window_sl = [15 50];

Tmax = 10.5;
Ttrunc = 0.5;

bouts = importdata(fullfile(paths_out.dataset, 'bouts.mat'));
bouts = bouts_formatting(bouts, thresholds);

paths_out = path_generator('folder', '/spontaneous_process', 'bouts_id', id_code, 'imfirst', false);
bouts_le = data_parser_new(bouts, 'period', 'loom', 'window', 'le', 'type', 'immobility', 'nloom', 2:20, 'min_dur', Ttrunc * 60);

T = bouts_le;

% % --- Dynamically find the subject-level group name from your configuration ---
% subject_level_struct_idx = find(cellfun(@(s) isfield(s, 'is_subject_level') && s.is_subject_level, random_effects_structure));
% if isempty(subject_level_struct_idx)
%     error('Plotting requires one grouping variable to be tagged with ''is_subject_level = true''.');
% end
% subject_level_group_name = random_effects_structure{subject_level_struct_idx}.group;

% =========================================================================
% 2. EXTRACT DATA FROM THE TABLE
% =========================================================================


% --- Get the names of the fields ---
re_coeff_names = {'nloom', 'avg_fs_1s', 'moving_flies', 'sloom', 'avg_sm', 'avg_ss', ...
    'time_since_last', 'n_generated_freezes', 'cum_freeze_time', 'avg_history_dur'};
matrix = bouts_le{:, re_coeff_names};
[rows, columns] = find(isnan(matrix));
unique(rows);
unique(columns);
matrix(unique(rows), :) = [];

% --- Get the Covariance Matrix ---
% This is the raw covariance matrix for this random effect group.
% cov_matrix = cov(matrix);

% --- Convert Covariance to Correlation ---
% A correlation matrix is much easier to interpret than a covariance matrix
% because all values are standardized to be between -1 and 1.
[corr_matrix] = corrcoef(matrix);

% =========================================================================
% 3. CREATE THE HEATMAP PLOT
% =========================================================================

figure('color', 'w', 'Name', 'Random Effects Correlation Matrix');
re_coeff_names = cellfun(@(x) strrep(x, '_', '\_'), re_coeff_names, 'UniformOutput', false);
re_coeff_names = cellfun(@(x) strrep(x, '\_\_', '*'), re_coeff_names, 'UniformOutput', false);

h = heatmap(re_coeff_names, re_coeff_names, corr_matrix);

% --- Format the plot for clarity ---
h.Title = sprintf('Random Effects Correlation Matrix for Group: ''%s''', subject_level_group_name);

% Setting color limits from -1 to 1 is crucial for interpreting correlations.
h.ColorLimits = [-1, 1];

% using 'redbluecmap' from the Mapping Toolbox (diverging colormap).
h.Colormap = redbluecmap;

h.XLabel = 'Random Effect';
h.YLabel = 'Random Effect';

% Adjust text size if labels overlap
if length(re_coeff_names) > 8
    h.FontSize = 8;
end

disp('Correlation matrix plot created successfully.');
