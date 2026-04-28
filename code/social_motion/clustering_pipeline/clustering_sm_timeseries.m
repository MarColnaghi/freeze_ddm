
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
bouts_proc = data_parser_new(bouts, 'type', 'immobility', 'period', 'loom', 'window', 'le', 'nloom', 2:20, 'min_dur', 30);

%% Extract the pre-freezing social motion
sm_freeze = extract_sm_from_bouts(bouts_proc, 'type', 'onlyfreeze', 'output_type', 'mat');
inizio = 0;
fine = 630;
sm_freeze_full = extract_sm_from_bouts(bouts_proc, 'type', 'onsets', 'output_type', 'mat', 'window', [inizio fine]);

% sm_freeze_full is trials x timesteps (each row is a trial, each column a
% timestep in the duration)

clustering_type = 'max';
paths = path_generator('folder', fullfile('social_motion','clustering', clustering_type), 'bouts_id', id_code, 'imfirst', false);

n_clusters = 10;

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
            threshold_sm = 0.6;
            [max_vals, metric] = max(centered_sm_freeze, [], 2);
            metric = double(metric);

            failed_to_reach = max_vals < threshold_sm;
            trunc_point = 30;
            early_peaks = metric < trunc_point;

        else
            % Metric: Center of Mass
            pow = 6;
            t = 1:size(centered_sm_freeze, 2);
            metric = sum(t .* (centered_sm_freeze.^pow), 2) ./ sum(centered_sm_freeze.^pow, 2);

            failed_to_reach = false(size(metric));
        end

        % --- Split valid (with a max between truncation and censoring)
        valid_mask = ~failed_to_reach & ~early_peaks;

        % No Peak
        invalid_mask = failed_to_reach;

        % Early Peak
        early_mask = early_peaks;

        % Metric is the index of the peak, valid metric is the same but
        % just for the valid trials

        valid_metric = metric(valid_mask);
        valid_indices = find(valid_mask);

        [~, sort_valid] = sort(valid_metric);
        valid_indices_sorted = valid_indices(sort_valid);

        n_valid = numel(valid_indices_sorted);

        n_clusters_eff = min(n_clusters, n_valid);

        edges = round(linspace(1, n_valid+1, n_clusters_eff+1));

        temp_cluster_ids = zeros(n_valid,1);

        for i = 1:n_clusters_eff
            temp_cluster_ids(edges(i):edges(i+1)-1) = i;
        end

        idx = zeros(size(metric));

        % early peaks → cluster 1
        idx(early_peaks) = 1;

        % valid clusters shift by +1
        idx(valid_indices_sorted) = temp_cluster_ids + 1;

        % failed → last cluster
        idx(failed_to_reach) = n_clusters_eff + 2;

        actual_clusters = n_clusters_eff + 2;

        % --- FINAL SORT (cluster first, then metric)
        metric(invalid_mask) = inf; % ensures failed go last
        [~, sort_order] = sortrows([idx, metric], [1 2]);

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

% Sort the matrices
sorted_matrix = centered_sm_freeze(sort_order, :);
bouts_sorted = bouts_proc(sort_order, :);

% Update the representation loop to account for the extra cluster
all_unique_clusters = unique(idx);
all_unique_clusters(all_unique_clusters == 0) = []; % Clean up
repr = nan(length(all_unique_clusters), size(centered_sm_freeze, 2));

for i = 1:length(all_unique_clusters)
    curr_id = all_unique_clusters(i);
    rows_in_cluster = (idx == curr_id);
    repr(i, :) = mean(sm_freeze_full(rows_in_cluster, :), 1);
end

% Define boundaries for plotting (using the sorted cluster IDs)

cluster_id_sorted = idx(sort_order);
boundaries = [0; find(diff(cluster_id_sorted)) + 0.5; size(centered_sm_freeze, 1)];

col = cmapper([], length(all_unique_clusters));
if any(invalid_mask)
    col.pca(length(all_unique_clusters),:) = [0.25 0.25 0.25];
end

% Visualize the Heatmap
fh = figure('color','w','Position',[100, 100, 600, 4000]);
tiledlayout(4, 1, 'Padding', 'compact', 'TileSpacing', 'tight')
nexttile(2, [3 1])
ax_distr(2) = gca;

hold on
x_axis = inizio:fine;
y_axis = 1:size(sorted_matrix, 1);
h = imagesc(x_axis, y_axis, sm_freeze_full(sort_order, :), [0 1.2]);
% set(h, 'AlphaData', ~isnan(sm_freeze(sort_order, :)));

scatter(bouts_proc.durations(sort_order), 1:height(bouts_proc), 2, '|', 'k')
apply_generic(gca, 'xlim', [0 630], 'no_y', true, 'ylim', [-100 size(centered_sm_freeze, 1) + 100], 'xtick', 0:120:600)
colormap(cbrewer2('Reds',[]));

%cb = colorbar(ax, 'Location', 'southoutside', 'FontSize', 18, 'LineWidth', 2);
%cb.Label.String = 'Social Motion';

for idx_cluster = 1:length(all_unique_clusters)
    fill([ax_distr(2).XLim(1)-2,ax_distr(2).XLim(1)-2,ax_distr(2).XLim(1)-17,ax_distr(2).XLim(1)-17], [boundaries(idx_cluster) boundaries(idx_cluster + 1) boundaries(idx_cluster + 1)   boundaries(idx_cluster)], ...
        col.pca(idx_cluster,:), 'EdgeColor','none', 'Clipping', 'off');
    text(mean([ax_distr(2).XLim(1)-2,ax_distr(2).XLim(1)-17]), mean([boundaries(idx_cluster); boundaries(idx_cluster + 1)]), num2str(idx_cluster),...
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
ax_distr(1) = gca;
hold on
for idx_cluster = 1:length(all_unique_clusters)
    plot(inizio:fine, repr(idx_cluster, :), 'LineWidth', 1.5, 'Color', col.pca(idx_cluster, :));
end
apply_generic(gca, 'ylim', [0 1.2],'ytick', [0 1.2], 'xlim', [0 630], 'no_x', true, 'font_size', 18);
ylabel({'Mean', 'Social Motion'})

linkaxes([ax_distr(:)], 'x')

exporter(fh, paths, 'clusters_profiles.pdf')
exporter(fh, paths, 'clusters_profiles.png')

axes(ax_distr(2))
clim([90 91.2])
exporter(fh, paths, 'clusters_profiles_nocolor.pdf')
exporter(fh, paths, 'clusters_profiles_nocolor.png')

axes(ax_distr(2))
clim([0 1.2])
% Similaritiesc

fh = figure('color', 'w','Position', [100, 100, 400, 400]);
D = pdist(repr, 'correlation');
Z = linkage(D, 'ward');
dh = dendrogram(Z, 'Reorder', 1:length(all_unique_clusters));
set(dh, 'Color', 'k', 'LineWidth', 2);
apply_generic(gca)
exporter(fh, paths, 'clusters_similarity.pdf')

fh_distr = figure('color', 'w','Position', [100, 100, 900, 550]);
tiledlayout(ceil(length(all_unique_clusters)/4), 4, 'Padding', 'loose', 'TileSpacing', 'compact')

fps = 60;
bin_size =  10 / fps;
t_vec = (inizio:fine) ./ fps;

for idx_cluster = 1:length(all_unique_clusters)
    ax_distr(idx_cluster) = nexttile;
    hold(ax_distr(idx_cluster), 'on')

    plot(t_vec, repr(idx_cluster, :), 'Color', col.pca(idx_cluster,:), 'LineWidth', 2)
    [m_med, i_med] = max(repr(idx_cluster, :));
    t_med = t_vec(i_med);
    text(t_med, m_med, num2str(round(t_med, 2)), 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center')

    % Histogram
    censor_thresh = 10.5;
    bouts_proc.durations_s(bouts_proc.durations_s > 10.5) = 10.75;

    edges = [min(bouts_proc.durations_s):bin_size:censor_thresh, 11];

    hh = histogram(ax_distr(idx_cluster), ...
        bouts_proc.durations_s(idx == idx_cluster), ...
        edges, ...
        'Normalization', 'pdf', ...
        'FaceColor', col.pca(idx_cluster,:), ...
        'EdgeColor', 'k', ...
        'DisplayStyle', 'bar', ...
        'LineWidth', 0.8);

    apply_generic(ax_distr(idx_cluster), 'xlim', [0 10], 'ylim', [-0.2 2.2], 'tick_length', 0.035, 'xticks', 0:2:10)
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
    xlim(ax_distr(idx_cluster), [0 11])
    ylim(ax_distr(idx_cluster), [-0.15 1.65])

    % --- Axis styling ---
    ax_distr(idx_cluster).Clipping = 'on';
    set(ax_distr(idx_cluster), 'Layer', 'top')
    ax_distr(idx_cluster).YAxis.Visible = 'off';

    % --- 5. Move xlabel dynamically ---
    ax_distr(idx_cluster).XLabel.Position(2) = xlabel_y_pos;

    %     % --- Censored data (>10.5 s) ---
    %     censor_thresh = 10.5;
    %     durations = bouts_proc.durations_s(idx == idx_cluster);
    %
    %     n_censored = sum(durations > censor_thresh);
    %
    %     if n_censored > 0
    %         % Convert count to "pdf-like" height (same normalization as histogram)
    %         bin_width = bin_size;
    %         total_n = numel(durations);
    %         cens_height = n_censored / (total_n * bin_width);
    %
    %         % Plot as a bar at the edge (slightly beyond 10)
    %         x_cens = 10.2;
    %
    %         bar(ax_distr(idx_cluster), x_cens, cens_height, ...
    %             'FaceColor', col.pca(idx_cluster,:), ...
    %             'EdgeColor', 'k', ...
    %             'BarWidth', 0.3);
    %     end

end

exporter(fh_distr, paths, 'clusters_freezedurations.pdf')
exporter(fh_distr, paths, 'clusters_freezedurations.png')

fh = figure('color','w','Position',[100, 100, 400, 400]);
hold on

for idx_cluster =  1:length(all_unique_clusters)

    [f, x] = ecdf(bouts_proc.durations_s(idx == idx_cluster));

    plot(x, f, 'Color', col.pca(idx_cluster, :), 'LineWidth', 2.5)

end

apply_generic(gca, 'xlim', [0 10], 'ylim', [0 1], 'tick_length', 0.025, 'xticks', 0:2:10)
ylabel('ecdf')
xlabel('Duration (s)')
exporter(fh, paths, 'ecdf.pdf')
exporter(fh, paths, 'ecdf.png')

%%
fh = figure('color', 'w','Position', [100, 100, 1100, 550]);
tlo = tiledlayout(ceil(n_clusters/4), 4, 'Padding', 'loose', 'TileSpacing', 'compact');

for idx_cluster = 1:n_clusters
    indices_of_cluster = find(idx == idx_cluster);
    peak_buffer = 30;

    % Initialize your data stacks
    data_stacks = {[], [], []}; % Before, During, After

    for i = 1:length(indices_of_cluster)
        row_idx = indices_of_cluster(i);
        current_signal = sm_freeze_full(row_idx, :);
        bout_end_idx = bouts_proc.durations(row_idx);
        [~, peak_time_idx] = max(current_signal);

        if bout_end_idx < (peak_time_idx - peak_buffer)
            data_stacks{1} = [data_stacks{1}; current_signal];
        elseif bout_end_idx > (peak_time_idx + peak_buffer)
            data_stacks{3} = [data_stacks{3}; current_signal];
        else
            data_stacks{2} = [data_stacks{2}; current_signal];
        end
    end

    nexttile
    hold on;
    colors = cbrewer2('Set1', 3);
    colors = colors([2, 1, 3], :);
    labels = {'Before', 'During', 'After'};

    scatter(600, 1.5, 500, col.pca(idx_cluster, :), 'filled')
    text(600, 1.5, num2str(idx_cluster), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'FontSize', 18)

    for j = 1:3
        current_group = data_stacks{j};
        if size(current_group, 1) > 1
            avg = mean(current_group, 1);
            sem = std(current_group, 0, 1) ./ sqrt(size(current_group, 1));

            % --- USE YOUR FUNCTION HERE ---
            % x, y1 (upper), y2 (lower), where (empty for all), varargin
            fill_between(x_axis, avg + sem, avg - sem, [], ...
                'FaceColor', colors(j,:), 'EdgeColor', 'none', 'FaceAlpha', 0.25);

            % Plot the mean line on top
            plot(x_axis, avg, 'Color', colors(j,:), 'LineWidth', 1.5);

        elseif size(current_group, 1) == 1
            plot(x_axis, current_group, 'Color', colors(j,:), 'LineWidth', 1.5);
        end
    end

    apply_generic(gca, 'ylim', [0 2], 'xlim', [min(x_axis) max(x_axis)], 'xticks', 0:120:600, 'no_y', true)
    xticklabels(0:2:10)
end

exporter(fh, paths, 'clusters_before_during_after.pdf')
exporter(fh, paths, 'clusters_before_during_after.png')


%% Check profiles vs Control Variables
columns = 6;
fh = figure('color','w','Position',[100, 100, 1600, 4000]);
tiledlayout(4, columns, 'Padding', 'compact', 'TileSpacing', 'tight')
nexttile(columns + 1, [3 3])
ax_distr(2) = gca;

hold on
x_axis = inizio:fine;
y_axis = 1:size(sorted_matrix, 1);
h = imagesc(x_axis, y_axis, sm_freeze_full(sort_order, :), [0 2.2]);
% set(h, 'AlphaData', ~isnan(sm_freeze(sort_order, :)));

% 2. Calculate Individual Peaks (Adjusted for x_axis start)
% Find the index of the max in each row
[~, peak_idx_raw] = max(sm_freeze_full(sort_order, :), [], 2);
% Shift by 'inizio' to align with the plot's X-coordinates
individual_peaks = x_axis(1) + (peak_idx_raw - 1);

% 1. Define your boundaries (ensure they are column vectors)
left_bound  = (individual_peaks - peak_buffer)';
right_bound = (individual_peaks + peak_buffer)';
y_axis_vec  = y_axis(:)'; % Ensure row vector for fliplr

% 2. Before Area Polygon
% From 'inizio' to 'left_bound'
x_before = [ones(size(left_bound)) * inizio, fliplr(left_bound)];
y_before = [y_axis_vec, fliplr(y_axis_vec)];

% 3. After Area Polygon
% From 'right_bound' to 'fine'
x_after = [right_bound, ones(size(right_bound)) * fine];
y_after = [y_axis_vec, fliplr(y_axis_vec)];

% 4. Plot the Shaded Areas
hold on;
% Before Area (using Red from Set1 as per your category colors)
fill(x_before, y_before, colors(1,:), ...
    'FaceAlpha', 0.3, 'EdgeColor', 'none', 'DisplayName', 'Before Zone');

% After Area (using Green from Set1 as per your category colors)
fill(x_after, y_after, colors(3,:), ...
    'FaceAlpha', 0.3, 'EdgeColor', 'none', 'DisplayName', 'After Zone');

% plot(individual_peaks - peak_buffer, y_axis,'LineWidth', 1, 'Color', [colors(1,:)]);
% plot(individual_peaks + peak_buffer, y_axis, 'LineWidth', 1, 'Color', [colors(3,:)]);
scatter(bouts_proc.durations(sort_order), 1:height(bouts_proc), 2, '|', 'k')


apply_generic(gca, 'xlim', [0 630], 'no_y', true, 'ylim', [-100 size(centered_sm_freeze, 1) + 100], 'xtick', 0:120:600)
colormap(cbrewer2('Reds',[]));

%cb = colorbar(ax, 'Location', 'southoutside', 'FontSize', 18, 'LineWidth', 2);
%cb.Label.String = 'Social Motion';

for idx_cluster = 1:length(all_unique_clusters)
    fill([ax_distr(2).XLim(1)-2,ax_distr(2).XLim(1)-2,ax_distr(2).XLim(1)-17,ax_distr(2).XLim(1)-17], [boundaries(idx_cluster) boundaries(idx_cluster + 1) boundaries(idx_cluster + 1)   boundaries(idx_cluster)], ...
        col.pca(idx_cluster,:), 'EdgeColor','none', 'Clipping', 'off');
    text(mean([ax_distr(2).XLim(1)-2,ax_distr(2).XLim(1)-17]), mean([boundaries(idx_cluster); boundaries(idx_cluster + 1)]), num2str(idx_cluster),...
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
ax_distr(1) = gca;
hold on
for idx_cluster = 1:length(all_unique_clusters)
    plot(inizio:fine, repr(idx_cluster, :), 'LineWidth', 1.5, 'Color', col.pca(idx_cluster, :));
end
apply_generic(gca, 'ylim', [0 1.2],'ytick', [0 1.2], 'xlim', [0 630], 'no_x', true, 'font_size', 18);
ylabel({'Mean', 'Social Motion'})

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
    values = zeros(length(all_unique_clusters), n_bins);
    for c = 1:length(all_unique_clusters)
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
    colors2 = cbrewer2(cmap_name, n_bins);
    for k = 1:numel(bh)
        bh(k).FaceColor = 'flat';
        bh(k).CData = colors2(k, :);
    end

    % 5. Alignment & Formatting
    ax_bar = gca;
    set(ax_bar, 'YLim', [0.5 length(all_unique_clusters) + 0.5]);
    apply_generic(ax_bar, 'no_x', true, 'yticks', 1:12);
end

exporter(fh, paths, 'clusters_profiles_with_controlvars.pdf')
exporter(fh, paths, 'clusters_profiles_with_controlvars.png')

%% --- 1. VALID MASK


valid_sorted = valid_mask(sort_order);   % already in sorted order

% --- 2. -------- PEAK ALIGNMENT --------

ref_peak = 315;
shifts_peak = round(ref_peak - individual_peaks);


% --- 3. -------- LINEAR FIT ALIGNMENT --------

x = individual_peaks(valid_sorted);
y = find(valid_sorted);

coeffs = polyfit(x, y, 1);
slope = coeffs(1);
intercept = coeffs(2);

n_rows = size(sorted_matrix,1);
n_cols = size(sorted_matrix,2);

y_all = (1:n_rows)';
x_fit_all = (y_all - intercept) / slope;

ref_fit = mean(x_fit_all(valid_sorted));
shifts_fit = round(ref_fit - x_fit_all);

% --- 4. -------- PAD MATRICES (NO CROPPING) --------

max_shift = ceil(max(abs([shifts_peak(:); shifts_fit(:)]), [], 'omitnan'));
pad = max_shift + 30;

new_cols = n_cols + 2*pad;

aligned_peak = nan(n_rows, new_cols);
aligned_fit  = nan(n_rows, new_cols);

base_idx = pad + (1:n_cols);

for i = 1:n_rows
    if valid_sorted(i)
        if ~isnan(shifts_peak(i))
            aligned_peak(i, base_idx + shifts_peak(i)) = sorted_matrix(i,:);
        end

        if ~isnan(shifts_fit(i))
            aligned_fit(i, base_idx + shifts_fit(i)) = sorted_matrix(i,:);
        end
    end
end


% --- 5. -------- CENTERED TIME AXIS --------

ref_col_peak = pad + round(ref_peak);
ref_col_fit  = pad + round(ref_fit);

x_axis_peak = (1:new_cols) - ref_col_peak;
x_axis_fit  = (1:new_cols) - ref_col_fit;


% --- 6. -------- BREAK TIMES (CENTERED) --------

break_times = bouts_proc.durations(sort_order);

break_times_peak = break_times + shifts_peak - round(ref_peak);
break_times_fit  = break_times + shifts_fit  - round(ref_fit);


% --- 7. -------- FILTER VALID TRIALS --------

sorted_matrix_valid = sorted_matrix(valid_sorted, :);
aligned_peak_valid  = aligned_peak(valid_sorted, :);
aligned_fit_valid   = aligned_fit(valid_sorted, :);

break_times_valid       = break_times(valid_sorted);
break_times_peak_valid  = break_times_peak(valid_sorted);
break_times_fit_valid   = break_times_fit(valid_sorted);

n_rows_valid = sum(valid_sorted);


% --- 10. -------- SORT BY PEAK MAGNITUDE --------

peak_mag = max(aligned_peak_valid, [], 2, 'omitnan');
[peak_mag_sorted, sort_idx_mag] = sort(peak_mag, 'descend');

aligned_peak_mag = aligned_peak_valid(sort_idx_mag, :);
break_times_peak_mag = break_times_peak_valid(sort_idx_mag);

n_rows_mag = size(aligned_peak_mag,1);

% --- Z-SCORE (row-wise)
mu = mean(aligned_peak_mag, 2, 'omitnan');
sigma = std(aligned_peak_mag, 0, 2, 'omitnan');
sigma(sigma == 0) = NaN;

z_peak = (aligned_peak_mag - mu) ./ sigma;

n_mag_clusters = 6;
N = numel(peak_mag_sorted);
group_size = ceil(N / n_mag_clusters);

idx_mag_cluster_sorted = zeros(N,1);

for k = 1:n_mag_clusters
    inds = (k-1)*group_size + 1 : min(k*group_size, N);
    idx_mag_cluster_sorted(inds) = k;
end

edges = find(diff(idx_mag_cluster_sorted) ~= 0);
boundaries_sm = [0; edges; N];

%% --- 11. -------- PLOTTING --------

new_boundaries = boundaries - boundaries(2);
new_boundaries(1) = [];
new_boundaries(end) = [];

fh = figure('color','w','Position',[100 100 1800 800]);
tiledlayout(4, 3, 'TileSpacing', 'compact', 'Padding', 'compact')
colsm = cbrewer2('Reds', 1);

bin_size = 4;

% ================= TOP ROW =================

% --- ORIGINAL
nexttile
ax_distr(1) = gca;
hold on
histogram(break_times_valid, 0:bin_size:90000, 'Normalization','pdf', ...
    'FaceColor','k','EdgeColor','none')

apply_generic(ax_distr(1), 'ylim', [0 .02], 'xlim', [0 630], 'font_size', 18);
ylabel('Break Density')

yyaxis right
avg = mean(sorted_matrix_valid,1,'omitnan');
plot(1:n_cols, avg, 'r-', 'LineWidth', 2, 'Color', colsm)
xline(1,'k--')
ylim([0 3])
set(gca,'XTick',[], 'FontSize', 18)
ylabel('Mean SM')
ax_distr(1).YAxis(2).Color = [0.78, 0.21, 0.24];

% --- PEAK ALIGNED
nexttile
ax_distr(2) = gca;
hold on
histogram(break_times_peak_valid, -601:bin_size:max(break_times_peak_valid), ...
    'Normalization','pdf', 'FaceColor','k','EdgeColor','none')

apply_generic(ax_distr(2), 'ylim', [0 .01], 'xlim', [-630 630]./2, 'font_size', 18);

yyaxis right

for idx_valid_cluster = 1:n_clusters
    avg = mean(aligned_peak_valid(temp_cluster_ids == idx_valid_cluster,:), 1,'omitnan');
    plot(x_axis_peak, avg, 'r-', 'LineWidth', 1.5, 'Color', col.pca(idx_valid_cluster + 1 , :))
end

avg = mean(aligned_peak_valid,1,'omitnan');
plot(x_axis_peak, avg, 'k-', 'LineWidth', 2, 'Color', 'k')

xline(1,'k--')
ylim([0 3])
set(gca,'XTick',[], 'FontSize', 18)
ax_distr(2).YAxis(2).Color = [0.78, 0.21, 0.24];

% --- MAGNITUDE SORTED
nexttile
ax_distr(3) = gca;
hold on
colmag = flipud(cbrewer2('Set2', n_mag_clusters));

for k = 1:n_mag_clusters

    data = break_times_peak_mag(idx_mag_cluster_sorted == k);
    data = data(~isnan(data));

    [f, x] = ecdf(data);

    plot(x, f, ...
        'Color', colmag(k,:), ...
        'LineWidth', 2);
end

ylabel('CDF')
set(gca,'FontSize',18)
apply_generic(ax_distr(3), 'ylim', [0 1], 'xlim', [-630 630]./2, 'font_size', 18);

yyaxis right
avg = mean(aligned_peak_mag,1,'omitnan');
for idx_mag_clusters = 1:n_mag_clusters
    plot(x_axis_peak, mean(aligned_peak_mag(idx_mag_cluster_sorted == idx_mag_clusters, :), 1, 'omitnan'), 'k-', 'Color', colmag(idx_mag_clusters,:), 'LineWidth', 1)
end

xline(1,'k--')
ylim([0 3])
set(gca,'XTick',[], 'FontSize', 18)
ax_distr(3).YAxis(2).Color = colsm;


% ================= HEATMAPS =================

% --- ORIGINAL
nexttile(4,[3 1])
ax_distr(4) = gca;
imagesc(1:n_cols, 1:n_rows_valid, sorted_matrix_valid, [0 3.2])
set(gca,'YDir','normal')
hold on
scatter(break_times_valid, 1:n_rows_valid, 2, '|', 'k')
apply_generic(gca, 'no_y', true, 'xlim', [0 630])
colormap(cbrewer2('Reds',[]))

for idx_cluster = 1:n_clusters
    fill([ax_distr(4).XLim(1)-2,ax_distr(4).XLim(1)-2,ax_distr(4).XLim(1)-17,ax_distr(4).XLim(1)-17], [new_boundaries(idx_cluster)  new_boundaries(idx_cluster +1) new_boundaries(idx_cluster +1)   new_boundaries(idx_cluster)], ...
        col.pca(idx_cluster,:), 'EdgeColor','none', 'Clipping', 'off');
    text(mean([ax_distr(4).XLim(1)-2,ax_distr(4).XLim(1)-17]), mean([new_boundaries(idx_cluster); new_boundaries(idx_cluster + 1)]), num2str(idx_cluster + 1),...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
    hold on

end

xlabel('Time');

% Add horizontal lines to separate clusters
hold on;
for idx_cluster = 1:n_clusters
    yline(new_boundaries(idx_cluster), 'k--', 'LineWidth', 1.1);
end


% --- PEAK ALIGNED
nexttile(5,[3 1])
ax_distr(5) = gca;
imagesc(x_axis_peak, 1:n_rows_valid, aligned_peak_valid, [0 3.2])
set(gca,'YDir','normal')
hold on
scatter(break_times_peak_valid, 1:n_rows_valid, 2, '|', 'k')

apply_generic(gca, 'no_y', true, 'xlim', [-630 630]./2)
colormap(cbrewer2('Reds',[]))

for idx_cluster = 1:n_clusters
    fill([ax_distr(5).XLim(1)-2,ax_distr(5).XLim(1)-2,ax_distr(5).XLim(1)-17,ax_distr(5).XLim(1)-17], [new_boundaries(idx_cluster)  new_boundaries(idx_cluster+1) new_boundaries(idx_cluster+1)   new_boundaries(idx_cluster)], ...
        col.pca(idx_cluster,:), 'EdgeColor','none', 'Clipping', 'off');
    text(mean([ax_distr(5).XLim(1)-2,ax_distr(5).XLim(1)-17]), mean([new_boundaries(idx_cluster ); new_boundaries(idx_cluster +1)]), num2str(idx_cluster +1 ),...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
    hold on

end

xlabel('Time');

% Add horizontal lines to separate clusters
hold on;
for idx_cluster = 1:n_clusters
    yline(new_boundaries(idx_cluster), 'k--', 'LineWidth', 1.1);
end

% --- MAG SORTED
nexttile(6,[3 1])
ax_distr(6) = gca;
imagesc(x_axis_peak, 1:n_rows_mag, aligned_peak_mag, [0 3.2])
set(gca,'YDir','normal')
hold on
scatter(break_times_peak_mag, 1:n_rows_mag, 2, '|', 'k')
apply_generic(gca, 'no_y', true, 'xlim', [-630 630]./2)
colormap(cbrewer2('Reds',[]))

x_left = ax_distr(6).XLim(1);

for k = 1:n_mag_clusters

    y1 = boundaries_sm(k) + 1;
    y2 = boundaries_sm(k + 1);

    fill([x_left-2, x_left-2, x_left-17, x_left-17], ...
        [y1 y2 y2 y1], ...
        colmag(k,:), ...
        'EdgeColor','none', ...
        'Clipping','off');

    text(x_left-9.5, mean([y1 y2]), num2str(k), ...
        'HorizontalAlignment','center', ...
        'VerticalAlignment','middle');
end

xlabel('Time');

% ================= LINK =================

linkaxes([ax_distr(1) ax_distr(4)], 'x')
linkaxes([ax_distr(2) ax_distr(5)], 'x')
linkaxes([ax_distr(3) ax_distr(6)], 'x')

exporter(fh, paths, 'peak_vs_magnitude_sorted.pdf')
exporter(fh, paths, 'peak_vs_magnitude_sorted.png')

% window = 15;
% axes(ax(2))


%% Now we plot the social motion timeseries based on whether the freeze ended or not at the peak
fh = figure('color','w','Position',[100 100 750 400]);
tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact')
nexttile
hold on

before = break_times_peak_valid < -window;
after  = break_times_peak_valid > 15;
during = ~(after | before);

% Pack conditions into a cell array
conds = {before, during, after};

for i = 1:3

    idxc = conds{i};
    data = aligned_peak_valid(idxc, :);

    % Mean and SEM
    mu  = mean(data, 1, 'omitnan');
    sem = std(data, 0, 1, 'omitnan') ./ sqrt(sum(~isnan(data), 1));

    % Valid mask
    m = ~isnan(sem);

    % Shaded error ONLY for "before" (to match your original code)
    if i == 1
        fill_between(x_axis_peak(m), ...
            mu(m) + sem(m), ...
            mu(m) - sem(m), ...
            [], ...
            'FaceColor', colors(i, :), ...
            'EdgeColor', 'none', ...
            'FaceAlpha', 0.25);
    end

    % Mean line
    plot(x_axis_peak, mu, ...
        'Color', colors(i, :), ...
        'LineWidth', 2);
end

apply_generic(gca, 'xlim', [-360 120]./2, 'xticks', -180:60:180, 'ylim', [0 2]);
xticklabels(-3:1:3)

ylabel('Social Motion')
xlabel('Time (Aligned to Peak)')

nexttile
hold on

sm_freeze = extract_sm_from_bouts(bouts_proc, 'type', 'onlyfreeze', 'output_type', 'mat', 'align', 'offset');
sm_freeze_valid = sm_freeze(sort_order(valid_sorted), :);
durs = bouts_proc.durations(sort_order(valid_sorted));

before_early = before & durs >= 30 & durs < 60;
before_mid = before & durs >= 60 & durs < 180;
before_late = before & durs >= 180;

sm_freeze_valid_beforeearly = sm_freeze_valid(before_early,:);
sm_freeze_valid_beforemid = sm_freeze_valid(before_mid,:);
sm_freeze_valid_beforelate = sm_freeze_valid(before_late,:);
sm_freeze_valid_during = sm_freeze_valid(during,:);
sm_freeze_valid_after = sm_freeze_valid(after,:);

sm_sets = {
    sm_freeze_valid_beforeearly
    sm_freeze_valid_beforemid
    sm_freeze_valid_beforelate
    sm_freeze_valid_during
    sm_freeze_valid_after
    };

colors_before = cbrewer2('Blues', 4);
colors_loop = [colors_before(2:end,:); colors(2:3, :)];


for i = 1:5

    data = sm_sets{i};

    % Mean and SEM
    mu = mean(data, 1, 'omitnan');
    sem = std(data, 0, 1, 'omitnan') ./ sqrt(sum(~isnan(data), 1));

    % Valid indices
    m = ~isnan(sem);

    % Time vector (aligned like your original intent)
    t = -numel(mu) + 1 : 0;

    % Shaded error
    fill_between(t(m), ...
        mu(m) + sem(m), ...
        mu(m) - sem(m), ...
        [], ...
        'FaceColor', colors_loop(i, :), ...
        'EdgeColor', 'none', ...
        'FaceAlpha', 0.25);

    % Mean line
    plot(t, mu, ...
        'Color', colors_loop(i, :), ...
        'LineWidth', 2);
    %     plot(t, sum(~isnan(data), 1), 'Color',  colors_loop(i, :))

end

zoomY = [0.3 .8];
zoomX = [-120 0];

fill([zoomX fliplr(zoomX)], [zoomY(1) zoomY(1) zoomY(2) zoomY(2)], [0 0 0],'FaceColor', 'none', 'EdgeColor', 'k', 'LineWidth', 2)
apply_generic(gca, 'xlim', [-300 0], 'xticks', -360:120:0, 'ylim', [0 2])
xticklabels(-6:2:0)
xlabel('Time (Aligned to Offset)')

% Get position of second tile
ax2 = gca;
pos = ax2.Position;

% Define inset relative to that tile
inset_pos = [
    pos(1) + 0.2*pos(3), ...
    pos(2) + 0.55*pos(4), ...
    0.5*pos(3), ...
    0.35*pos(4)
    ];

ax_inset = axes('Position', inset_pos);
hold(ax_inset, 'on')
box(ax_inset, 'on')

for i = 1:3

    data = sm_sets{i};
    if isempty(data)
        continue
    end

    mu  = mean(data, 1, 'omitnan');
    sem = std(data, 0, 1, 'omitnan') ./ sqrt(sum(~isnan(data), 1));

    m = ~isnan(sem);
    t = -numel(mu) + 1 : 0;

    fill_between(t(m), ...
        mu(m) + sem(m), ...
        mu(m) - sem(m), ...
        [], ...
        'FaceColor', colors_loop(i,:), ...
        'EdgeColor', 'none', ...
        'FaceAlpha', 0.25);

    plot(t, mu, ...
        'Color', colors_loop(i,:), ...
        'LineWidth', 1.5);

    apply_generic(ax_inset, 'xlim', zoomX, 'xticks', -360:120:0, 'ylim', zoomY, 'yticks', [0.3 .8])
    xticklabels(-6:2:0)

end

exporter(fh, paths, 'sm_aligned_2_offset.pdf')
exporter(fh, paths, 'sm_aligned_2_offset.png')

%%

fh = figure('Position', [100 100 1300 700], 'Color', 'w');
tl = tiledlayout(3, 4, 'TileSpacing', 'loose', 'Padding', 'compact');
i = 0;
bin_size = 0.1;


for idx_sm = 1:3
    for idx_ls = 0:1
        for idx_fs = 1:2

            nexttile
            hold on
            i = i + 1;

            [freezes_quant, mask_quant] = quantilizer_v2(bouts_proc, 'indexed_quantile', struct('sm', idx_sm, 'ls', idx_ls, 'fs', idx_fs));

            sum_freezes(i) = numel(mask_quant);
            histogram(freezes_quant.durations_s, min(freezes_quant.durations_s):bin_size:max(freezes_quant.durations_s), 'Normalization', 'pdf');

            ax_main = gca;

            apply_generic(ax_main, 'ylim', [0 2], 'xlim', [0 10.5], 'no_y', true, 'font_size', 18)

            % --- Create inset inside current tile ---

            pos = ax_main.Position;

            inset_pos = [
                pos(1) + 0.2*pos(3), ...
                pos(2) + 0.80*pos(4), ...
                0.40*pos(3), ...
                0.30*pos(4)
                ];

            ax_inset = axes('Position', inset_pos);
            hold(ax_inset, 'on')
            box(ax_inset, 'on')


            % Define bins (same bin_size or coarser for clarity)
            edges = 0.5:1:max(all_unique_clusters) + 0.5;  % coarser bins work better in small inset
            counts = histcounts(idx(mask_quant), edges) ./ numel(mask_quant);

            b = bar(ax_inset, counts, 'FaceColor', 'flat', 'EdgeColor', 'none');

            % Assign one color per bar
            b.CData = col.pca(1:numel(counts), :);
            yline(0.2, 'k--', 'LineWidth', 0.5)

            % --- Formatting inset ---
            xlim(ax_inset, [0 13])   % zoom on short durations
            ylim(ax_inset, [0 .35])

            apply_generic(ax_inset, 'xlim', [0 13], 'ylim', [0 .35], 'font_size', 14, 'no_x', true, 'line_width', 1)
            hold(ax_inset, 'off')

            axes(ax_main)
            data = sm_freeze(mask_quant, :);

            mu  = mean(data, 1, 'omitnan');
            sem = std(data, 0, 1, 'omitnan') ./ sqrt(sum(~isnan(data), 1));

            t_full = (1:numel(mu)) ./ fps;
            m = ~isnan(sem);

            yyaxis right
            % --- SHADED SEM ---
            fill_between(t_full(m), ...
                mu(m) + sem(m), ...
                mu(m) - sem(m), ...
                [], ...
                'FaceColor', [0.8 0.24 0.24], ... % or use a color per condition
                'EdgeColor', 'none', ...
                'FaceAlpha', 0.2);

            % --- MEAN TRACE ---
            plot(t_full, mu, 'Color', [0.8 0.24 0.24], 'LineWidth', 1.5)

            set(gca, 'FontSize', 18)
            ylim([0 2])
            ylabel('Social Motion')
            xlabel('Time (s)')

        end
    end
end


exporter(fh, paths, 'sm_divided_by_conditions.pdf')
exporter(fh, paths, 'sm_divided_by_conditions.png')


%%

% =========================
% --- USER PARAMETERS -----
% =========================

group_var_name = 'ls';     % 'ls' OR 'fs'
group_mode     = 'raw';    % 'raw' OR 'median'

% =========================
% --- FIGURE SETUP --------
% =========================

fh_distr = figure('Color','w','Position',[100 100 900 550]);

n_clusters = numel(all_unique_clusters);
n_cols = 4;
n_rows = ceil(n_clusters / n_cols);

tiledlayout(n_rows, n_cols, 'Padding','loose','TileSpacing','compact');

% =========================
% --- CONSTANTS -----------
% =========================

fps = 60;
bin_size = 18 / fps;
t_vec = (inizio:fine) ./ fps;

censor_thresh = 10.5;
edges = min(bouts_proc.durations_s):bin_size:censor_thresh;
edges(end+1) = 11;

% --- Preprocess durations ONCE ---
durations = bouts_proc.durations_s;
durations_capped = durations;
durations_capped(durations_capped > censor_thresh) = 10.75;

% =========================
% --- GROUPING LOGIC ------
% =========================

group_var = bouts_proc.(group_var_name);

switch group_mode

    case 'raw'
        unique_groups = unique(group_var(~isnan(group_var)));
        group_masks = arrayfun(@(g) group_var == g, unique_groups, 'UniformOutput', false);
        group_labels = arrayfun(@(g) sprintf('%s = %.2f', group_var_name, g), ...
                                unique_groups, 'UniformOutput', false);

    case 'median'
        med = median(group_var, 'omitnan');

        group_masks = {
            group_var <= med
            group_var >  med
        };

        group_labels = {
            sprintf('%s <= med', group_var_name)
            sprintf('%s > med',  group_var_name)
        };

    otherwise
        error('Unknown group_mode')
end

n_groups = numel(group_masks);

% --- Colors ---
coli = cmapper([], 5);
cols = coli.vars.(group_var_name);
cols = cols([4, 6], :);

% =========================
% --- MAIN LOOP -----------
% =========================

for idx_cluster = 1:n_clusters

    ax = nexttile;
    hold(ax, 'on')

    mask_cluster = (idx == idx_cluster);

    % -------------------------
    % LOOP OVER GROUPS
    % -------------------------
    for i_g = 1:n_groups

        maski = mask_cluster & group_masks{i_g};

        if sum(maski) < 5
            continue
        end

        % --- Mean trace ---
        mean_trace = mean(sm_freeze_full(maski, :), 1);

        plot(ax, t_vec, mean_trace, ...
            'Color', cols(i_g,:), ...
            'LineWidth', 1.3);

        % --- Histogram ---
        histogram(ax, durations_capped(maski), edges, ...
            'Normalization','pdf', ...
            'DisplayStyle','stairs', ...
            'EdgeColor', cols(i_g,:), ...
            'LineWidth', 1.25);
    end

    % -------------------------
    % REPRESENTATIVE PEAK
    % -------------------------
    [m_med, i_med] = max(repr(idx_cluster, :));
    t_med = t_vec(i_med);

    text(ax, t_med, m_med, sprintf('%.2f', t_med), ...
        'VerticalAlignment','bottom', ...
        'HorizontalAlignment','center');

    % -------------------------
    % HEATMAP STRIP
    % -------------------------
    imagesc(ax, t_vec, -0.52, repr(idx_cluster,:), [0 1.2]);
    colormap(ax, cbrewer2('Reds', []));

    ax.YAxis.Visible = 'off';
    ax.Clipping = 'on';
    set(ax, 'Layer', 'top');

    % -------------------------
    % X LABEL LOGIC
    % -------------------------
    if idx_cluster == 9
        xlabel(ax, 'Freeze Duration (s)', 'FontSize', 14);
    else
        ax.XTickLabel = [];
    end

    apply_generic(ax, 'xlim', [0 11], 'ylim', [-0.2 1.5], 'tick_length', 0.035, 'xticks', 0:2:10, 'font_size', 18) 
    xtickangle(0)
end

% =========================
% --- EXPORT --------------
% =========================

exporter(fh_distr, paths, ...
    sprintf('clusters_freezedurations_%s_%s.pdf', group_var_name, group_mode));

exporter(fh_distr, paths, ...
    sprintf('clusters_freezedurations_%s_%s.png', group_var_name, group_mode));

%%

fh = figure('color','w','Position',[100 100 1800 800]);
tiledlayout(4, 3, 'TileSpacing', 'compact', 'Padding', 'compact')
group_var_name = 'ls';     % 'ls' OR 'fs'

bin_size = 4;

% ================= LOOM SPEED SETUP =================
unique_ls = unique(bouts_proc.ls);
n_ls = numel(unique_ls);
coli = cmapper([], 5);
cols = coli.vars.(group_var_name);
cols = cols([4, 6], :);

% ================= TOP ROW =================

% --- ORIGINAL
nexttile
clear ax_distr
ax_distr(1) = gca;
hold on

histogram(break_times_valid, 0:bin_size:90000, ...
    'Normalization','pdf', ...
    'FaceColor','k','EdgeColor','none')

apply_generic(ax_distr(1), 'ylim', [0 .02], 'xlim', [0 630], 'font_size', 18);
ylabel('Break Density')

yyaxis right
avg = mean(sorted_matrix_valid,1,'omitnan');
n_cols = size(sorted_matrix,2);

plot(1:n_cols, avg, 'LineWidth', 2, 'Color', [0.78, 0.21, 0.24 0.4])

xline(1,'k--')
ylim([0 3])
set(gca,'XTick',[], 'FontSize', 18)
ylabel('Mean SM')
ax_distr(1).YAxis(2).Color = [0.78, 0.21, 0.24];

% --- PEAK ALIGNED (NOW SPLIT BY LOOM SPEED)
nexttile
ax_distr(2) = gca;
hold on

edges_peak = -601:bin_size:max(break_times_peak_valid);
apply_generic(ax_distr(2), 'ylim', [0 .02], 'xlim', [-630 630]./2, 'font_size', 18);

for i_ls = 1:n_ls

    ls_val = unique_ls(i_ls);
    mask = bouts_proc.ls(sort_order(valid_mask)) == ls_val;

    % --- Histogram ---
    histogram(break_times_peak_valid(mask), ...
        edges_peak, ...
        'Normalization','pdf', ...
        'DisplayStyle','stairs', ...
        'EdgeColor', cols(i_ls,:), ...
        'LineWidth', 1.5);



end
apply_generic(ax_distr(2), 'ylim', [0 .02], 'xlim', [-630 630]./2, 'font_size', 18);

yyaxis right


for i_ls = 1:n_ls

    % --- Mean trace ---
    avg = mean(aligned_peak_valid(mask, :), 1, 'omitnan');

    plot(x_axis_peak, avg, ...
        'LineWidth', 2, ...
        'Color', [cols(i_ls,:)]); % transparent

end

yyaxis right
xline(0,'k--')   % peak-aligned → 0 is correct
ylim([0 3])

set(gca,'XTick',[], 'FontSize', 18)
ax_distr(2).YAxis(2).Color = [0.78, 0.21, 0.24];

% Optional legend
legend(arrayfun(@(x) sprintf('LS = %d', x), unique_ls, 'UniformOutput', false), ...
    'Location','northeast')

% --- MAGNITUDE SORTED
nexttile
ax_distr(3) = gca;
hold on

colsm = flipud(cbrewer2('Set2', n_mag_clusters));

for k = 1:n_mag_clusters
    data = break_times_peak_mag(idx_mag_cluster_sorted == k);
    data = data(~isnan(data));

    [f, x] = ecdf(data);

    plot(x, f, ...
        'Color', colsm(k,:), ...
        'LineWidth', 2);
end

ylabel('CDF')
set(gca,'FontSize',18)
apply_generic(ax_distr(3), 'ylim', [0 1], 'xlim', [-630 630]./2, 'font_size', 18);

yyaxis right
for idx_mag_clusters = 1:n_mag_clusters
    plot(x_axis_peak, ...
        mean(aligned_peak_mag(idx_mag_cluster_sorted == idx_mag_clusters, :), 1, 'omitnan'), 'k-',...
        'Color', colsm(idx_mag_clusters,:), ...
        'LineWidth', .5)
end

xline(1,'k--')
ylim([0 3])
set(gca,'XTick',[], 'FontSize', 18)
ax_distr(3).YAxis(2).Color = [0.78, 0.21, 0.24];

% ================= HEATMAPS =================

% --- ORIGINAL
nexttile(4,[3 1])
ax_distr(4) = gca;

imagesc(1:n_cols, 1:n_rows_valid, sorted_matrix_valid, [0 3.2])
set(gca,'YDir','normal')
hold on

scatter(break_times_valid, 1:n_rows_valid, 2, '|', 'k')

apply_generic(gca, 'no_y', true, 'xlim', [0 630])
colormap(cbrewer2('Reds',[]))

for idx_cluster = 1:n_clusters
    fill([ax_distr(4).XLim(1)-2,ax_distr(4).XLim(1)-2,ax_distr(4).XLim(1)-17,ax_distr(4).XLim(1)-17], ...
        [new_boundaries(idx_cluster) new_boundaries(idx_cluster+1) new_boundaries(idx_cluster+1) new_boundaries(idx_cluster)], ...
        col.pca(idx_cluster,:), 'EdgeColor','none', 'Clipping','off');

    text(mean([ax_distr(4).XLim(1)-2,ax_distr(4).XLim(1)-17]), ...
        mean([new_boundaries(idx_cluster) new_boundaries(idx_cluster+1)]), ...
        num2str(idx_cluster+1), ...
        'HorizontalAlignment','center', 'VerticalAlignment','middle');
end

xlabel('Time');

for idx_cluster = 1:n_clusters
    yline(new_boundaries(idx_cluster), 'k--', 'LineWidth', 1.1);
end

% --- PEAK ALIGNED (UNCHANGED HEATMAP)
nexttile(5,[3 1])
ax_distr(5) = gca;

imagesc(x_axis_peak, 1:n_rows_valid, aligned_peak_valid, [0 3.2])
set(gca,'YDir','normal')
hold on

scatter(break_times_peak_valid, 1:n_rows_valid, 2, '|', 'k')

apply_generic(gca, 'no_y', true, 'xlim', [-630 630]./2)
colormap(cbrewer2('Reds',[]))

for idx_cluster = 1:n_clusters
    fill([ax_distr(5).XLim(1)-2,ax_distr(5).XLim(1)-2,ax_distr(5).XLim(1)-17,ax_distr(5).XLim(1)-17], ...
        [new_boundaries(idx_cluster) new_boundaries(idx_cluster+1) new_boundaries(idx_cluster+1) new_boundaries(idx_cluster)], ...
        col.pca(idx_cluster,:), 'EdgeColor','none', 'Clipping','off');

    text(mean([ax_distr(5).XLim(1)-2,ax_distr(5).XLim(1)-17]), ...
        mean([new_boundaries(idx_cluster) new_boundaries(idx_cluster+1)]), ...
        num2str(idx_cluster+1), ...
        'HorizontalAlignment','center', 'VerticalAlignment','middle');
end

xlabel('Time');

for idx_cluster = 1:n_clusters
    yline(new_boundaries(idx_cluster), 'k--', 'LineWidth', 1.1);
end

% --- MAG SORTED
nexttile(6,[3 1])
ax_distr(6) = gca;

imagesc(x_axis_peak, 1:n_rows_mag, aligned_peak_mag, [0 3.2])
set(gca,'YDir','normal')
hold on

scatter(break_times_peak_mag, 1:n_rows_mag, 2, '|', 'k')

apply_generic(gca, 'no_y', true, 'xlim', [-630 630]./2)
colormap(cbrewer2('Reds',[]))

x_left = ax_distr(6).XLim(1);

for k = 1:n_mag_clusters
    y1 = boundaries_sm(k) + 1;
    y2 = boundaries_sm(k + 1);

    fill([x_left-2, x_left-2, x_left-17, x_left-17], ...
        [y1 y2 y2 y1], ...
        colsm(k,:), ...
        'EdgeColor','none', ...
        'Clipping','off');

    text(x_left-9.5, mean([y1 y2]), num2str(k), ...
        'HorizontalAlignment','center', ...
        'VerticalAlignment','middle');
end

xlabel('Time');

% ================= LINK =================

linkaxes([ax_distr(1) ax_distr(4)], 'x')
linkaxes([ax_distr(2) ax_distr(5)], 'x')
linkaxes([ax_distr(3) ax_distr(6)], 'x')

%exporter(fh, paths, 'peak_vs_magnitude_sorted_ls_split.pdf')
%exporter(fh, paths, 'peak_vs_magnitude_sorted_ls_split.png')
% %%
% fh = figure('color','w','Position',[100, 100, 400, 400]);
% hold on
%
% for idx_cluster =  1:length(all_unique_clusters)
%
%     [f, x] = ecdf(bouts_proc.durations_s(idx == idx_cluster));
%
%     f_tvec = interp1(x(2:end), f(2:end), t_vec);
%     plot3(t_vec, f_tvec,  repr(idx_cluster,:),  'Color', col.pca(idx_cluster, :), 'LineWidth', 2.5);
%
% end
%
% apply_generic(gca, 'xlim', [0 10.5], 'ylim', [0 1], 'tick_length', 0.025, 'xticks', 0:2:10)
% ylabel('ecdf')
% xlabel('Duration (s)')
%
% % Set up the video writer
% v = VideoWriter(fullfile(paths.fig, 'ecdf3d.avi'));
% open(v);
% % Specify the number of frames for the animation
% numFrames = 400;
%
% % Rotate the surface plot and capture frames
% for k = 1:numFrames
%
%     if k < 60
%         view([0, 90]);
%     else
%
%         % Change the view angle
%         view([k, 90]);
%     end
%     % Capture the current frame
%     frame = getframe(gcf);
%     writeVideo(v, frame);
% end
% % Complete the video writing process
% close(v);
% % Close the figure display
% close(gcf);
%
% % Here we fit the different clusters