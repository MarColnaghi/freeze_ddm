clearvars

% Load the table first. We will take advantage of an already existing
% dataset.
threshold_imm = 2; threshold_mob = 2; threshold_pc = 4; id_code = sprintf('imm%d_mob%d_pc%d', threshold_imm, threshold_mob, threshold_pc);
paths = path_generator('folder', 'social_motion/clustering', 'bouts_id', id_code, 'imfirst', false);
bouts = importdata(fullfile(paths.dataset, 'bouts.mat'));
motion_cache = importdata(fullfile(paths.cache_path, 'motion_cache.mat'));

% Here you can change the window of loom-evoked
thresholds = define_thresholds('le_window', struct('le_window_sl', [5 55], 'le_window_fl', [5 55]));
thresholds = define_thresholds;

bouts = bouts_formatting(bouts, thresholds);

% Process the data: set thresholds for sm, fs and ln. Set minimum duration.
bouts_proc = data_parser_new(bouts, 'type', 'immobility', 'period', 'loom', 'window', 'le', 'nloom', 2:20, 'min_dur', 30, 'max_dur', 630);

%% Extract the pre-freezing social motion
sm_freeze = extract_sm_from_bouts(bouts_proc, 'type', 'onlyfreeze', 'output_type', 'mat');
inizio = 0;
fine = 630;
sm_freeze_full = extract_sm_from_bouts(bouts_proc, 'type', 'onsets', 'output_type', 'mat', 'window', [inizio fine]);

% sm_freeze_full is trials x timesteps (each row is a trial, each column a
% timestep in the duration)

% centered_sm_freeze = sm_freeze_full - mean(sm_freeze_full, 2);
% centered_sm_freeze = zscore(sm_freeze_full, [], 2);

clustering_type = 'max';
paths = path_generator('folder', fullfile('social_motion','clustering', clustering_type), 'bouts_id', id_code, 'imfirst', false);

n_clusters = 8;
col = cmapper([], n_clusters);

switch clustering_type

    case {'max'}
        centered_sm_freeze = sm_freeze;
        if strcmp(clustering_type, 'max')
            % Metric: Time of Max Peak
            [~, metric] = max(centered_sm_freeze, [], 2);
            total_freeze_dur = bouts_proc.durations;
            relative_peaks = metric./total_freeze_dur;
        else
            % Metric: Center of Mass
            pow = 6;
            t = 1:size(centered_sm_freeze, 2);
            metric = sum(t .* (centered_sm_freeze.^ pow), 2) ./ sum(centered_sm_freeze.^ pow, 2);
        end
        
        % 1. Get the sort order based on the metric
        [~, sort_order] = sort(relative_peaks);
        
        % 2. Assign Clusters based on position in the sorted list
        n = size(centered_sm_freeze, 1);
        % Create cluster labels (e.g., 1,1,1, 2,2,2...)
        temp_cluster_ids = ceil((1:n)' * n_clusters / n);
        
        % 3. Map these cluster labels back to the original row indices
        idx = zeros(n, 1);
        idx(sort_order) = temp_cluster_ids;

    case 'first_cross'
        threshold = 0.8;
        centered_sm_freeze = sm_freeze;

        metric = nan(size(centered_sm_freeze,1),1);
        for row = 1:size(centered_sm_freeze,1)
            tmp = find(centered_sm_freeze(row,:) > threshold, 1, 'first');
            if ~isempty(tmp)
                metric(row) = tmp;
            end
        end

        % 1. Get the sort order based on the metric
        [~, sort_order] = sort(metric);
        
        % 2. Assign Clusters based on position in the sorted list
        n = size(centered_sm_freeze, 1);
        % Create cluster labels (e.g., 1,1,1, 2,2,2...)
        temp_cluster_ids = ceil((1:n)' * n_clusters / n);
        
        % 3. Map these cluster labels back to the original row indices
        idx = zeros(n, 1);
        idx(sort_order) = temp_cluster_ids;

end


% Define boundaries for plotting (using the sorted cluster IDs)
cluster_id_sorted = idx(sort_order); 
boundaries = [0; find(diff(cluster_id_sorted)) + 0.5; size(centered_sm_freeze, 1)];

% Sort the matrices
sorted_matrix = centered_sm_freeze(sort_order, :);
bouts_sorted = bouts_proc(sort_order, :);

% Extract mean timeseries
repr = nan(n_clusters, size(centered_sm_freeze, 2));
for idx_cluster = 1:n_clusters
    % This now works because idx always refers to original row indices
    rows_in_cluster = (idx == idx_cluster);
    if any(rows_in_cluster)
        repr(idx_cluster, :) = mean(sm_freeze(rows_in_cluster, :), 1);
    end
end

% Visualize the Heatmap
fh = figure('color','w','Position',[100, 100, 600, 4000]);
tiledlayout(4, 1, 'Padding', 'compact', 'TileSpacing', 'tight')
nexttile(1, [4 1])
ax(2) = gca;

hold on
x_axis = (inizio:fine) + 1; 
y_axis = 1:size(sorted_matrix, 1);
h = imagesc(x_axis, y_axis, sm_freeze(sort_order, :), [0 3.0]);
% set(h, 'AlphaData', ~isnan(sm_freeze(sort_order, :)));

scatter(bouts_proc.durations(sort_order), 1:height(bouts_proc), 2, '|', 'k')
apply_generic(gca, 'xlim', [0 630], 'no_y', true, 'ylim', [-100 size(centered_sm_freeze, 1) + 100], 'xtick', 0:120:600)
colormap(cbrewer2('Reds',[]));

%cb = colorbar(ax, 'Location', 'southoutside', 'FontSize', 18, 'LineWidth', 2);
%cb.Label.String = 'Social Motion';

for idx_cluster = 1:n_clusters
    fill([ax(2).XLim(1)-2,ax(2).XLim(1)-2,ax(2).XLim(1)-17,ax(2).XLim(1)-17], [boundaries(idx_cluster) boundaries(idx_cluster + 1) boundaries(idx_cluster + 1)   boundaries(idx_cluster)], ...
        col.pca(idx_cluster,:), 'EdgeColor','none', 'Clipping', 'off');
    text(mean([ax(2).XLim(1)-2,ax(2).XLim(1)-17]), mean([boundaries(idx_cluster); boundaries(idx_cluster + 1)]), num2str(idx_cluster),...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
    hold on

end

xlabel('Time');

% Add horizontal lines to separate clusters
hold on;
for idx_cluster = 1:length(boundaries)
    yline(boundaries(idx_cluster), 'k--', 'LineWidth', 1.1);
end


%%

% Visualize the Heatmap
fh = figure('color','w','Position',[100, 100, 600, 4000]);
tiledlayout(4, 1, 'Padding', 'compact', 'TileSpacing', 'tight')
nexttile(1, [4 1])
ax(2) = gca;

% Define the number of points for the normalized timeline
num_points = 300; 
query_points = linspace(0, 1, num_points);

% Pre-allocate the matrix for warped data
sm_freeze_warped = zeros(size(sorted_matrix, 1), num_points);

for idx_bout = 1:size(sorted_matrix, 1)
    % Get the actual duration of the current bout
    current_dur = bouts_proc.durations(idx_bout);
    
    % Create the original time vector for this specific bout (0 to 1)
    original_time = linspace(0, 1, current_dur);
    
    % Extract the motion data (assuming sm_freeze contains the raw vectors)
    % We take only the columns up to current_dur if it's padded with NaNs
    current_data = sm_freeze(idx_bout, 1:current_dur);
    
    % Interpolate the data onto the 100-point query grid
    sm_freeze_warped(idx_bout, :) = interp1(original_time, current_data, query_points, 'linear');
end

h = imagesc(query_points, y_axis, sm_freeze_warped(sort_order, :), [0 1.0]);

% set(h, 'AlphaData', ~isnan(sm_freeze(sort_order, :)));

%scatter(bouts_proc.durations(sort_order), 1:height(bouts_proc), 2, '|', 'k')
apply_generic(gca, 'xlim', [0 1], 'no_y', true, 'ylim', [-100 size(centered_sm_freeze, 1) + 100], 'xtick', 0:0.2:1)
colormap(cbrewer2('Reds',[]));

%cb = colorbar(ax, 'Location', 'southoutside', 'FontSize', 18, 'LineWidth', 2);
%cb.Label.String = 'Social Motion';

% for idx_cluster = 1:n_clusters
%     fill([ax(2).XLim(1)-2,ax(2).XLim(1)-2,ax(2).XLim(1)-17,ax(2).XLim(1)-17], [boundaries(idx_cluster) boundaries(idx_cluster + 1) boundaries(idx_cluster + 1)   boundaries(idx_cluster)], ...
%         col.pca(idx_cluster,:), 'EdgeColor','none', 'Clipping', 'off');
%     text(mean([ax(2).XLim(1)-2,ax(2).XLim(1)-17]), mean([boundaries(idx_cluster); boundaries(idx_cluster + 1)]), num2str(idx_cluster),...
%         'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
%     hold on
% 
% end

xlabel('Time');

% Add horizontal lines to separate clusters
hold on;
for idx_cluster = 1:length(boundaries)
    yline(boundaries(idx_cluster), 'k--', 'LineWidth', 1.1);
end