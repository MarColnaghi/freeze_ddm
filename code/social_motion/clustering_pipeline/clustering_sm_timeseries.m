
clearvars

% Load the table first. We will take advantage of an already existing
% dataset.
threshold_imm = 2; threshold_mob = 2; threshold_pc = 4; id_code = sprintf('imm%d_mob%d_pc%d', threshold_imm, threshold_mob, threshold_pc);
paths = path_generator('folder', 'social_motion/clustering', 'bouts_id', id_code, 'imfirst', false);
bouts = importdata(fullfile(paths.dataset, 'bouts.mat'));
motion_cache = importdata(fullfile(paths.cache_path, 'motion_cache.mat'));

% Here you can change the window of loom-evoked
thresholds = define_thresholds('le_window', struct('le_window_sl', [5 55], 'le_window_fl', [5 55]));
% thresholds = define_thresholds;

bouts = bouts_formatting(bouts, thresholds);

% Process the data: set thresholds for sm, fs and ln. Set minimum duration.
bouts_proc = data_parser_new(bouts, 'type', 'immobility', 'period', 'loom', 'window', 'le', 'nloom', 2:20, 'min_dur', 12);

%% Extract the pre-freezing social motion
sm_freeze = extract_sm_from_bouts(bouts_proc, 'type', 'onlyfreeze', 'output_type', 'mat');
inizio = 0;
fine = 60;
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
    case 'pca'
        % 1. Run PCA
        centered_sm_freeze = zscore(sm_freeze_full, [], 2);
        [~, score] = pca(centered_sm_freeze);
        
        % 2. Cluster
        pcs_to_use = 30;
        data_for_clustering = score(:, 1:pcs_to_use);
        [idx, ~, ~, diffs] = kmeans(data_for_clustering, n_clusters, 'Replicates', 5, 'MaxIter', 10000, 'Distance', 'correlation');
        
       % 1. Calculate current cluster means to find peak timing
        temp_repr = nan(n_clusters, size(centered_sm_freeze, 2));
        for i = 1:n_clusters
            temp_repr(i, :) = mean(centered_sm_freeze(idx == i, :), 1);
        end
        
        % 2. Find the time-index of the max peak for each cluster
        [~, peak_times] = max(temp_repr, [], 2);
        
        % 3. Create a mapping from OLD cluster IDs to NEW (sorted) IDs
        [~, cluster_rank] = sort(peak_times); % Which old ID is earliest?
        new_idx_map = zeros(n_clusters, 1);
        for i = 1:n_clusters
            new_idx_map(cluster_rank(i)) = i; 
        end
        
        % 4. Reassign idx to the new "earlier-to-later" numbering
        idx = new_idx_map(idx);
        
        % 5. Now calculate distance to own centroid (using new idx) for sorting within clusters
        dist_to_assigned = arrayfun(@(i) diffs(i, cluster_rank(new_idx_map == idx(i))), (1:length(idx))');
        
        % 6. Final Sort Order: By new Cluster ID, then by distance
        [~, sort_order] = sortrows([idx, dist_to_assigned], [1 2]);

    case {'max', 'com'}
        centered_sm_freeze = sm_freeze_full;
        if strcmp(clustering_type, 'max')
            % Metric: Time of Max Peak
            % 1. Define your threshold
            threshold_sm = 0.3;

            % 2. Find the maximum value and its index for each row
            % centered_sm_freeze is assumed to be your [Observations x Time] matrix
            [max_vals, metric] = max(centered_sm_freeze, [], 2);

            % 3. Identify which rows failed to reach the threshold
            % This creates a logical mask (1 for fail, 0 for pass)
            failed_to_reach = max_vals < threshold_sm;

            % 4. Apply the threshold condition
            % We set the index to NaN for any row where the max wasn't high enough.
            % Note: metric must be converted to double to hold NaN values.
            metric = double(metric);
            metric(failed_to_reach) = NaN;

            % --- OPTIONAL: Clean up the results ---
            % If you want to know how many observations were rejected:
            num_rejected = sum(failed_to_reach);
            fprintf('Excluded %d rows that did not exceed threshold %0.2f\n', num_rejected, threshold_sm);

        else
            % Metric: Center of Mass
            pow = 6;
            t = 1:size(centered_sm_freeze, 2);
            metric = sum(t .* (centered_sm_freeze.^ pow), 2) ./ sum(centered_sm_freeze.^ pow, 2);
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

    case 'first_cross'
        threshold = 0.2;
        centered_sm_freeze = sm_freeze_full;

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
        repr(idx_cluster, :) = mean(sm_freeze_full(rows_in_cluster, :), 1);
    end
end

% Visualize the Heatmap
fh = figure('color','w','Position',[100, 100, 600, 4000]);
tiledlayout(4, 1, 'Padding', 'compact', 'TileSpacing', 'tight')
nexttile(2, [3 1])
ax(2) = gca;

hold on
x_axis = inizio:fine; 
y_axis = 1:size(sorted_matrix, 1);
h = imagesc(x_axis, y_axis, sm_freeze_full(sort_order, :), [0 2.2]);
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

nexttile(1)
ax(1) = gca;
hold on
for idx_cluster = 1:n_clusters
    plot(inizio:fine, repr(idx_cluster, :), 'LineWidth', 1.5, 'Color', col.pca(idx_cluster, :));
end
apply_generic(gca, 'ylim', [0 1.2],'ytick', [0 1.2], 'xlim', [0 630], 'no_x', true, 'font_size', 18);
ylabel({'Mean', 'Social Motion'})

linkaxes([ax(:)], 'x')

exporter(fh, paths, 'clusters_profiles.pdf')
exporter(fh, paths, 'clusters_profiles.png')

% axes(ax(2))
% clim([90 91.2])
% exporter(fh, paths, 'clusters_profiles_nocolor.pdf')
% exporter(fh, paths, 'clusters_profiles_nocolor.png')

% Similarities

fh = figure('color', 'w','Position', [100, 100, 400, 400]);
D = pdist(repr, 'correlation');
Z = linkage(D, 'ward');
dh = dendrogram(Z, 'Reorder', 1:n_clusters);
set(dh, 'Color', 'k', 'LineWidth', 2);
apply_generic(gca)
exporter(fh, paths, 'clusters_similarity.pdf')
 
fh_distr = figure('color', 'w','Position', [100, 100, 900, 550]);
tiledlayout(ceil(n_clusters/4), 4, 'Padding', 'loose', 'TileSpacing', 'compact')

fps = 60;
bin_size =  10 / fps;
t_vec = (inizio:fine) ./ fps;

for idx_cluster = 1:n_clusters
    ax_distr(idx_cluster) = nexttile;
    hold(ax_distr(idx_cluster), 'on')

    plot(t_vec, repr(idx_cluster, :), 'Color', col.pca(idx_cluster,:), 'LineWidth', 2)
    [m_med, i_med] = max(repr(idx_cluster, :));
    t_med = t_vec(i_med);
    text(t_med, m_med, num2str(round(t_med, 2)), 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center')  

    % Histogram
    histogram(ax_distr(idx_cluster), bouts_proc.durations_s(idx == idx_cluster), ...
        min(bouts_proc.durations_s):bin_size:max(bouts_proc.durations_s) + 1, ...
        'Normalization', 'pdf', ...
        'FaceColor', col.pca(idx_cluster,:), ...
        'EdgeColor', 'k', ...
        'DisplayStyle', 'bar',...
        'LineWidth', 0.8);

    apply_generic(ax_distr(idx_cluster), 'xlim', [0 10], 'ylim', [-0.2 2.2], 'tick_length', 0.025, 'xticks', 0:2:10)
    xtickangle(0)

    % Conditional xlabel (like your original logic)
    if idx_cluster == 9
        xlabel(ax_distr(idx_cluster), 'Freeze Duration (s)', 'FontSize', 18)
    else
        xlabel(ax_distr(idx_cluster), '')
        xticklabels(ax_distr(idx_cluster), [])
    end

    % --- 1. Y limits ---
    y_min = -0.5;
    y_max = 1;
    ylim(ax_distr(idx_cluster), [y_min, y_max]);
    y_range = y_max - y_min;

    % --- 2. Proportional layout ---
    img_y_pos     = y_min + (0.06 * y_range);
    xlabel_y_pos  = y_min - (0.15 * y_range);

    % --- 3. Prepare data for imagesc ---
    img_data = mean(sm_freeze_full(idx == idx_cluster, :));

    % --- 4. Plot imagesc with controlled position ---
    imagesc(ax_distr(idx_cluster), t_vec, -0.52, repr(idx_cluster,:), [0 1.2]);

    colormap(ax_distr(idx_cluster), cbrewer2('Reds', []))   % apply to axis only
    xlim(ax_distr(idx_cluster), [0 10])
    ylim(ax_distr(idx_cluster), [-0.15 1.65])

    % --- Axis styling ---
    ax_distr(idx_cluster).Clipping = 'on';
    set(ax_distr(idx_cluster), 'Layer', 'top')
    ax_distr(idx_cluster).YAxis.Visible = 'off';

    % --- 5. Move xlabel dynamically ---
    ax_distr(idx_cluster).XLabel.Position(2) = xlabel_y_pos;
end

exporter(fh_distr, paths, 'clusters_freezedurations.pdf')
exporter(fh_distr, paths, 'clusters_freezedurations.png')

fh = figure('color','w','Position',[100, 100, 400, 400]);
hold on

for idx_cluster =  1:n_clusters
    
    [f, x] = ecdf(bouts_proc.durations_s(idx == idx_cluster));

    plot(x, f, 'Color', col.pca(idx_cluster, :), 'LineWidth', 2.5)

end

apply_generic(gca, 'xlim', [0 10.5], 'ylim', [0 1], 'tick_length', 0.025, 'xticks', 0:2:10)
ylabel('ecdf')
xlabel('Duration (s)')
exporter(fh, paths, 'ecdf.pdf')
exporter(fh, paths, 'ecdf.png')

%% Check profiles vs Control Variables
columns = 6;
fh = figure('color','w','Position',[100, 100, 1600, 4000]);
tiledlayout(4, columns, 'Padding', 'compact', 'TileSpacing', 'tight')
nexttile(columns + 1, [3 3])
ax(2) = gca;

hold on
x_axis = inizio:fine; 
y_axis = 1:size(sorted_matrix, 1);
h = imagesc(x_axis, y_axis, sm_freeze_full(sort_order, :), [0 2.2]);
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

nexttile(1, [1 3])
ax(1) = gca;
hold on
for idx_cluster = 1:n_clusters
    plot(inizio:fine, repr(idx_cluster, :), 'LineWidth', 1.5, 'Color', col.pca(idx_cluster, :));
end
apply_generic(gca, 'ylim', [0 1.2],'ytick', [0 1.2], 'xlim', [0 630], 'no_x', true, 'font_size', 18);
ylabel({'Mean', 'Social Motion'})

linkaxes([ax(:)], 'x')

plot_configs = {
    bouts_proc.moving_flies, 'Greys',   5,  10;
    bouts_proc.sloom / 25,   'Blues',   2,  11;
    bouts_proc.nloom,        'Purples', 20, 12
};

for i = 1:size(plot_configs, 1)
    data_raw = plot_configs{i, 1};
    cmap_name = plot_configs{i, 2};
    n_bins    = plot_configs{i, 3};
    tile_idx  = plot_configs{i, 4};
    
    % 1. Pre-allocate and calculate histogram counts
    values = zeros(n_clusters, n_bins);
    for c = 1:n_clusters
        values(c, :) = histcounts(data_raw(idx == c), n_bins);
    end
    
    % 2. Normalize to Percentages (using element-wise division)
    perc_values = (values ./ sum(values, 2)) * 100;
    
    % 3. Plotting
    nexttile(tile_idx, [3 1])
    hold on
    xline([0 20 40 60 80 100], 'k--')
    bh = bar(perc_values, 'stacked', 'Horizontal', 'on', 'EdgeColor', 'black');
    
    % 4. Apply Colors
    colors = cbrewer2(cmap_name, n_bins);
    for k = 1:numel(bh)
        bh(k).FaceColor = 'flat';
        bh(k).CData = colors(k, :);
    end
    
    % 5. Alignment & Formatting
    ax_bar = gca;
    set(ax_bar, 'YLim', [0.5 n_clusters + 0.5]);
    apply_generic(ax_bar, 'no_x', true, 'yticks', 1:12);
end

exporter(fh, paths, 'clusters_profiles_with_controlvars.pdf')
exporter(fh, paths, 'clusters_profiles_with_controlvars.png')

%%
fh = figure('color','w','Position',[100, 100, 1600, 800]);

cluster = 5;
indices_of_cluster = find(idx == cluster);

for idx_bouts_in_cluster = 1:length(indices_of_cluster)
    
    hold on
    i = indices_of_cluster(idx_bouts_in_cluster);
    plot(sm_freeze_full(i, :) + idx_bouts_in_cluster, x_axis + 1, 'Color', [0.9 0.9 0.9]);
    plot(sm_freeze(i, :) + idx_bouts_in_cluster, 1:size(sm_freeze, 2), 'Color', 'k');

    [c, ts] = max(sm_freeze_full(i, :) + idx_bouts_in_cluster);
    plot(c, ts, 'Color', 'r', 'Marker', 'o');
end

apply_generic(gca, 'ylim', [0 630])



%%
fh = figure('color','w','Position',[100, 100, 400, 400]);
hold on

for idx_cluster =  1:n_clusters
    
    [f, x] = ecdf(bouts_proc.durations_s(idx == idx_cluster));

    f_tvec = interp1(x(2:end), f(2:end), t_vec);
    plot3(t_vec, f_tvec,  repr(idx_cluster,:),  'Color', col.pca(idx_cluster, :), 'LineWidth', 2.5);

end

apply_generic(gca, 'xlim', [0 10.5], 'ylim', [0 1], 'tick_length', 0.025, 'xticks', 0:2:10)
ylabel('ecdf')
xlabel('Duration (s)')

% Set up the video writer
v = VideoWriter(fullfile(paths.fig, 'ecdf3d.avi'));
open(v);
% Specify the number of frames for the animation
numFrames = 400;

% Rotate the surface plot and capture frames
for k = 1:numFrames

    if k < 60
        view([0, 90]);
    else

        % Change the view angle
        view([k, 90]);
    end
    % Capture the current frame
    frame = getframe(gcf);
    writeVideo(v, frame);
end
% Complete the video writing process
close(v);
% Close the figure display
close(gcf);

% Here we fit the different clusters