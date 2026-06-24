
clearvars

% Load the table first. We will take advantage of an already existing
% dataset.
threshold_imm = 2; threshold_mob = 2; threshold_pc = 4; id_code = sprintf('imm%d_mob%d_pc%d', threshold_imm, threshold_mob, threshold_pc);
paths = path_generator('folder', 'collision_analysis', 'bouts_id', id_code, 'imfirst', false);
bouts = importdata(fullfile(paths.dataset, 'bouts.mat'));
motion_cache = importdata(fullfile(paths.cache_path, 'motion_cache.mat'));

% Here you can change the window of loom-evoked
thresholds = define_thresholds('le_window', struct('le_window_sl', [5 55], 'le_window_fl', [5 55]));
% thresholds = define_thresholds;

bouts = bouts_formatting(bouts, thresholds);

% Process the data: set thresholds for sm, fs and ln. Set minimum duration.
bouts_proc = data_parser_new(bouts, 'type', 'immobility', 'period', 'loom', 'window', 'le', 'nloom', 2:20, 'min_dur', 30);
bouts_proc = impose_contact_threshold(bouts_proc, 'threshold', 80);

%% Extract the pre-freezing social motion
inizio = 0;
fine = 630;
distance_ili = extract_sm_from_bouts(bouts_proc, 'type', 'onsets', 'output_type', 'mat', 'window', [inizio fine], 'cache', 'mindist_cache');

% distance_ili is trials x timesteps (each row is a trial, each column a
% timestep in the duration)

clustering_type = 'first_cross';
paths = path_generator('folder', fullfile('collision_analysis', clustering_type), 'bouts_id', id_code, 'imfirst', false);

n_clusters = 10;

switch clustering_type
    case 'pca'
        % 1. Run PCA
        centered_sm_freeze = zscore(distance_ili, [], 2);
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

        centered_sm_freeze = distance_ili;

        if strcmp(clustering_type, 'max')
            threshold_sm = .5;
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
        threshold = 50;
        centered_sm_freeze = distance_ili;
        data = centered_sm_freeze;

        n = size(data,1);
        metric = nan(n,1);

        % --- 1. Compute metric (first crossing) ---
        for i = 1:n
            tmp = find(data(i,:) < threshold, 1, 'first');
            if ~isempty(tmp)
                metric(i) = tmp;
            end
        end

        failed_to_reach = isnan(metric);
        trunc_point = 30;
        early_peaks = metric < trunc_point;

        % --- 2. Define masks ---
        valid_mask  = ~failed_to_reach & ~early_peaks & ~isnan(metric);
        early_mask  = early_peaks;
        failed_mask = failed_to_reach;

        % --- 3. Sort valid trials ---
        valid_idx = find(valid_mask);
        [~, order_valid] = sort(metric(valid_mask));
        valid_sorted = valid_idx(order_valid);

        n_valid = numel(valid_sorted);

        % --- 4. Assign clusters to valid trials ---
        n_clusters_eff = min(n_clusters, n_valid);

        edges = round(linspace(1, n_valid+1, n_clusters_eff+1));

        cluster_valid = zeros(n_valid,1);
        for k = 1:n_clusters_eff
            cluster_valid(edges(k):edges(k+1)-1) = k;
        end

        % --- 5. Build final cluster index ---
        idx = nan(n,1);

        % early peaks → cluster 1
        idx(early_mask) = 1;

        % valid → clusters 2..(n_clusters_eff+1)
        idx(valid_sorted) = cluster_valid + 1;

        % failed → last cluster
        idx(failed_mask) = n_clusters_eff + 2;

        actual_clusters = n_clusters_eff + 2;

        % --- 6. Final sorting (cluster first, then metric) ---
        metric_for_sort = metric;
        metric_for_sort(isnan(metric_for_sort)) = inf;

        [~, sort_order] = sortrows([idx, metric_for_sort], [1 2]);
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
    repr(i, :) = mean(distance_ili(rows_in_cluster, :), 1);
end

% Define boundaries for plotting (using the sorted cluster IDs)

cluster_id_sorted = idx(sort_order);
boundaries = [0; find(diff(cluster_id_sorted)) + 0.5; size(centered_sm_freeze, 1)];

%%

col = cmapper([], length(all_unique_clusters));
if any(failed_mask)
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
h = imagesc(x_axis, y_axis, distance_ili(sort_order, :), [0 200]);
% set(h, 'AlphaData', ~isnan(sm_freeze(sort_order, :)));

scatter(bouts_proc.durations(sort_order), 1:height(bouts_proc), 2, '|', 'k')
apply_generic(gca, 'xlim', [0 630], 'no_y', true, 'ylim', [-100 size(centered_sm_freeze, 1) + 100], 'xtick', 0:120:600)
colormap(flipud(cbrewer2('Reds',[])));

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
apply_generic(gca, 'ylim', [0 300],'ytick', [0 60 120 300], 'xlim', [0 630], 'no_x', true, 'font_size', 18);
ylabel({'Mean', 'Min Distance'})
yline(60, 'k--')
linkaxes([ax_distr(:)], 'x')

exporter(fh, paths, 'clusters_profiles.pdf')
exporter(fh, paths, 'clusters_profiles.png')

axes(ax_distr(2))
clim([0 1])
exporter(fh, paths, 'clusters_profiles_nocolor.pdf')
exporter(fh, paths, 'clusters_profiles_nocolor.png')

axes(ax_distr(2))
clim([0 200])
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
    img_data = mean(distance_ili(idx == idx_cluster, :));

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

individual_peaks = x_axis(1) + (peak_idx_raw - 1);
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
break_times(break_times > 630) = NaN;

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
