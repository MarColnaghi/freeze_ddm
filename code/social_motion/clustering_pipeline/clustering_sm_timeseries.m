
clearvars

% Load the table first. We will take advantage of an already existing
% dataset.
threshold_imm = 2; threshold_mob = 2; threshold_pc = 4; id_code = sprintf('imm%d_mob%d_pc%d', threshold_imm, threshold_mob, threshold_pc);
paths = path_generator('folder', 'social_motion/clustering', 'bouts_id', id_code, 'imfirst', false);
bouts = importdata(fullfile(paths.dataset, 'bouts.mat'));
motion_cache = importdata(fullfile(paths.cache_path, 'motion_cache.mat'));

% Here you can change the window of loom-evoked
thresholds = define_thresholds('le_window', struct('le_window_sl', [5 55], 'le_window_fl', [5 55]));
%thresholds = define_thresholds;

bouts = bouts_formatting(bouts, thresholds);

% Process the data: set thresholds for sm, fs and ln. Set minimum duration.
bouts_proc = data_parser_new(bouts, 'type', 'immobility', 'period', 'loom', 'window', 'le', 'nloom', 2:20, 'min_dur', 30);

% Extract the pre-freezing social motion
sm_freeze = extract_sm_from_bouts(bouts_proc, 'type', 'onlyfreeze', 'output_type', 'mat');
sm_freeze_full = extract_sm_from_bouts(bouts_proc, 'type', 'onsets', 'output_type', 'mat', 'window', [0 630]);

%% sm_freeze_full is trials x timesteps (each row is a trial, each column a
% timestep in the duration)

% centered_sm_freeze = sm_freeze_full - mean(sm_freeze_full, 2);
% centered_sm_freeze = zscore(sm_freeze_full, [], 2);
centered_sm_freeze = sm_freeze_full;

% Sort by CoM
t = 1:size(centered_sm_freeze, 2);
pow = 6;
t_cm = sum(t .* (centered_sm_freeze.^ pow), 2) ./ sum(centered_sm_freeze.^ pow  , 2);
% t_cm = sum(t .* centered_sm_freeze, 2) ./ sum(centered_sm_freeze, 2);
[t_cm_sorted, t_cm_cluster_id] = sort(t_cm);

% Sort by Max
[max_sm, max_cluster_id] = max(centered_sm_freeze, [], 2);
[t_max_sorted, t_max_cluster_id] = sort(max_cluster_id);

% Visualize the Heatmap
fh = figure('color','w','Position',[100, 100, 600, 900]);
tiledlayout(4, 1, 'Padding', 'compact', 'TileSpacing', 'tight')
nexttile(2, [3 1])
hold on
sort_order = t_max_cluster_id;
h = imagesc(sm_freeze_full(sort_order, :), [0 1.2]);
set(h, 'AlphaData', ~isnan(sm_freeze(sort_order, :)));
scatter(bouts_proc.durations(sort_order), 1:height(bouts_proc), 3, '|', 'k')
ax = gca;
apply_generic(gca, 'xlim', [0 630], 'no_y', true, 'ylim', [-100 size(centered_sm_freeze, 1) + 100], 'xtick', 0:120:600)
colormap(cbrewer2('Reds',[]));
%cb = colorbar(ax, 'Location', 'southoutside', 'FontSize', 18, 'LineWidth', 2);
%cb.Label.String = 'Social Motion';

n_clusters = 10;

n = length(sort_order);
cluster_id = ceil((1:n)' * n_clusters / n);
boundaries = [0; find(diff(cluster_id)) + 0.5; size(centered_sm_freeze, 1)];
sorted_matrix = centered_sm_freeze(sort_order, :);
bouts_sorted = bouts_proc(sort_order, :);

col = cmapper([], n_clusters);

for idx_cluster = 1:n_clusters
    fill([ax.XLim(1)-2,ax.XLim(1)-2,ax.XLim(1)-17,ax.XLim(1)-17], [boundaries(idx_cluster) boundaries(idx_cluster + 1) boundaries(idx_cluster + 1)   boundaries(idx_cluster)], ...
        col.pca(idx_cluster,:), 'EdgeColor','none', 'Clipping', 'off');
    text(mean([ax.XLim(1)-2,ax.XLim(1)-17]), mean([boundaries(idx_cluster); boundaries(idx_cluster + 1)]), num2str(idx_cluster),...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
end

xlabel('Time');

% Add horizontal lines to separate clusters
hold on;
for idx_cluster = 1:length(boundaries)
    yline(boundaries(idx_cluster), 'k--', 'LineWidth', 1.1);
end

nexttile(1, [1 1])
repr = nan(n_clusters, size(centered_sm_freeze, 2));

for idx_cluster = 1:n_clusters
    hold on
    repr(idx_cluster, 1:size(centered_sm_freeze, 2)) = mean(sorted_matrix(cluster_id == idx_cluster, :));
    plot(1:size(centered_sm_freeze, 2), repr(idx_cluster, :), 'LineWidth', 1.5, 'Color', col.pca(idx_cluster, :));
end

ylabel('Median Social Motion')
apply_generic(gca, 'ylim', [0 1.2],'ytick', [0 1.2], 'no_x', true);
exporter(fh, paths, 'clusters_profiles.pdf')
exporter(fh, paths, 'clusters_profiles.png')

fh = figure('color','w','Position',[100, 100, 400, 400]);
D = pdist(repr, 'euclidean');
Z = linkage(D, 'ward');
dendrogram(Z);
exporter(fh, paths, 'clusters_similarity.pdf')
exporter(fh, paths, 'clusters_similarity.png')

fh = figure('color','w','Position',[100, 100, 900, 550]);
tiledlayout(3, 3, 'Padding', 'compact', 'TileSpacing', 'compact')

fps = 60;
bin_size =  10 / fps;
t_vec = 0:1/fps:10.5;

for idx_cluster = 1:n_clusters
    ax = nexttile;
    hold(ax, 'on')

    % Histogram
    histogram(ax, bouts_sorted.durations_s(cluster_id == idx_cluster), ...
        min(bouts_proc.durations_s):bin_size:900, ...
        'Normalization', 'pdf', ...
        'FaceColor', col.pca(idx_cluster,:), ...
        'EdgeColor', 'k', ...
        'DisplayStyle', 'bar',...
        'LineWidth', 1.5);

    apply_generic(ax, 'xlim', [0 10], 'ylim', [-0.2 2.2])

    plot(t_vec, repr(idx_cluster, :), 'Color', col.pca(idx_cluster,:), 'LineWidth', 2)


    % Conditional xlabel (like your original logic)
    if idx_cluster == 9
        xlabel(ax, 'Duration (s)')
    else
        xlabel(ax, '')
        xticklabels(ax, [])
    end

    % --- 1. Y limits ---
    y_min = -0.5;
    y_max = 1;
    ylim(ax, [y_min, y_max]);
    y_range = y_max - y_min;

    % --- 2. Proportional layout ---
    img_y_pos     = y_min + (0.06 * y_range);
    xlabel_y_pos  = y_min - (0.15 * y_range);

    % --- 3. Prepare data for imagesc ---
    img_data = mean(sm_freeze_full(cluster_id == idx_cluster, :));

    % --- 4. Plot imagesc with controlled position ---
    imagesc(ax, t_vec, -0.52, repr(idx_cluster,:), [0 1.2]);

    colormap(ax, cbrewer2('Reds', []))   % apply to axis only
    xlim(ax, [0 10])
    ylim(ax, [-0.15 1.65])

    % --- Axis styling ---
    ax.Clipping = 'on';
    set(ax, 'Layer', 'top')
    ax.YAxis.Visible = 'off';

    % --- 5. Move xlabel dynamically ---
    ax.XLabel.Position(2) = xlabel_y_pos;
end

exporter(fh, paths, 'clusters_freezedurations.pdf')
exporter(fh, paths, 'clusters_freezedurations.png')

fh = figure('color','w','Position',[100, 100, 400, 400]);
hold on

for idx_cluster =  1:n_clusters
    
    [f, x] = ecdf(bouts_sorted.durations_s(cluster_id == idx_cluster));

    plot(x, f, 'Color', col.pca(idx_cluster, :), 'LineWidth', 2.5)

end
apply_generic(gca, 'xlim', [0 10.5], 'ylim', [0 1])
ylabel('ecdf')
xlabel('Time')
exporter(fh, paths, 'ecdf.pdf')
exporter(fh, paths, 'ecdf.png')