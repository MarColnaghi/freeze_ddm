
clearvars

% Load the table first. We will take advantage of an already existing
% dataset.
threshold_imm = 2; threshold_mob = 2; threshold_pc = 4; id_code = sprintf('imm%d_mob%d_pc%d', threshold_imm, threshold_mob, threshold_pc);
paths = path_generator('folder', 'social_motion/pca', 'bouts_id', id_code, 'imfirst', false);
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
centered_sm_freeze = zscore(sm_freeze_full, [], 2);
% centered_sm_freeze = sm_freeze_full;

% centered_sm_freeze = sm_freeze_full;
[coeff, score, latent, tsquared, explained] = pca(centered_sm_freeze);

fh = figure('color','w','Position',[100, 100, 600, 650]);
tiledlayout(2, 1)

% 3. Visualize the Explained Variance (Scree Plot)
nexttile
pareto(explained);
cumvar = cumsum(explained);
pc_to_use = 1:find(cumvar > 90, 1);
xlabel('PCs');
ylabel('Variance Explained (%)');
apply_generic(gca)

% Plot the first three Principal Components over time
nexttile
plot(coeff(:, 1:5), 'LineWidth', 2);
legend('PC1','PC2','PC3', 'PC4','PC5', 'Location', 'eastoutside');
xlabel('Time (frames)');
ylabel('PC');
apply_generic(gca)
exporter(fh, paths, 'pcs.pdf')
exporter(fh, paths, 'pcs.png')

% Reconstruct using only the first 3 PCs
n_pc = 3;
reconstructed = score(:,1:n_pc) * coeff(:,1:n_pc)' + mean(centered_sm_freeze);
residuals = centered_sm_freeze - reconstructed;

fh = figure('color','w','Position',[100, 100, 600, 650]);
tiledlayout(2, 1)
nexttile
imagesc(residuals, [0 3]); % Look for structured patterns here
title(['Residuals using ' num2str(n_pc) ' PCs']);
apply_generic(gca)

nexttile
hold on

trial_id = 445;
n_pc = 30;
recon = score(:,1:n_pc) * coeff(:,1:n_pc)';
plot(centered_sm_freeze(trial_id,:)); hold on;
plot(recon(trial_id,:));
legend('Original','Reconstruction');
apply_generic(gca)

%% Kmeans clustering

% Determine how many PCs and clusters to have
% pc_to_use = 1:50;
% data_for_clustering = score(:, pc_to_use);
% 
% K_range = 1:30;
% within_sum = zeros(size(K_range));
% 
% for k = K_range
%     [~, ~, sumd] = kmeans(data_for_clustering, k, 'Replicates', 5, 'MaxIter',10000);
%     within_sum(k) = sum(sumd);
% end
% 
% figure;
% plot(K_range, within_sum, '-o');
% xlabel('Number of clusters');
% ylabel('Total within-cluster sum of squares');

n_clusters = 9;
pc_to_use = 1:30; cumvar(max(pc_to_use));
data_for_clustering = score(:, pc_to_use);
col = cmapper([], n_clusters);

[idx, C, sums, diffs] = kmeans(data_for_clustering, n_clusters, 'Replicates', 5, 'MaxIter', 10000, 'Distance', 'correlation');

repr = nan(n_clusters, size(sm_freeze_full, 2));

for idx_cluster = 1:n_clusters
    hold on
    repr(idx_cluster, 1:size(sm_freeze_full, 2)) = mean(sm_freeze_full(idx == idx_cluster, :));
    plot(1:size(sm_freeze_full, 2), repr(idx_cluster, :), 'LineWidth', 1.5, 'Color', col.pca(idx_cluster, :));
end

% Distance to the ASSIGNED cluster center for each trial
dist_to_assigned = zeros(size(idx));
for i = 1:length(idx)
    dist_to_assigned(i) = diffs(i, idx(i));
end

% 3. Visualize the clusters in PC space
fh = figure('color','w','Position',[100, 100, 600, 650]);
scatter3(score(:,1), score(:,2), score(:,3), 12, idx, 'filled');
xlabel('PC1'); ylabel('PC2'); zlabel('PC3');
grid on; colormap(col.pca);
apply_generic(gca)
exporter(fh, paths, 'clusters_in_pc_space.pdf')
exporter(fh, paths, 'clusters_in_pc_space.png')

% Multi-level sort: First by Cluster ID (idx), then by Distance
% We create a combined matrix to keep indices linked
combine_for_sort = [idx, dist_to_assigned, (1:length(idx))']; 
sorted_data = sortrows(combine_for_sort, [1 2]); 

% Extract the new order
sort_order = sorted_data(:, 3);
sorted_idx = sorted_data(:, 1);

boundaries = [0; find(diff(sorted_idx)) + 0.5; size(sm_freeze_full, 1)];

%% Visualize the Heatmap
fh = figure('color','w','Position',[100, 100, 600, 650]);
tiledlayout(3, 1, 'Padding', 'compact', 'TileSpacing', 'tight')
nexttile(2, [2 1])
hold on
imagesc(sm_freeze_full(sort_order, :), [0 1.2]);
scatter(bouts_proc.durations(sort_order), 1:height(bouts_proc), 3, '|', 'k')
ax = gca;
apply_generic(gca, 'xlim', [0 630], 'no_y', true, 'ylim', [-100 size(sm_freeze_full, 1) + 100], 'xtick', 0:120:600)
colormap(cbrewer2('Reds',[]));
cb = colorbar(ax, 'Location', 'southoutside', 'FontSize', 18, 'LineWidth', 2);
cb.Label.String = 'Social Motion';
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
repr = nan(n_clusters, size(sm_freeze_full, 2));

for idx_cluster = 1:n_clusters
    hold on
    repr(idx_cluster, 1:size(sm_freeze_full, 2)) = mean(sm_freeze_full(idx == idx_cluster, :));
    plot(1:size(sm_freeze_full, 2), repr(idx_cluster, :), 'LineWidth', 1.5, 'Color', col.pca(idx_cluster, :));
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

%%

fh = figure('color','w','Position',[100, 100, 900, 550]);
tiledlayout(3, 4, 'Padding', 'compact', 'TileSpacing', 'compact')

fps = 60;
bin_size =  10 / fps;
t_vec = 0:1/fps:10.5;

for idx_cluster = 1:n_clusters
    ax = nexttile;
    hold(ax, 'on')

    % Histogram
    histogram(ax, bouts_proc.durations_s(idx == idx_cluster), ...
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
    img_data = mean(sm_freeze_full(idx == idx_cluster, :));

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

%%

fh = figure('color','w','Position',[100, 100, 400, 400]);
hold on

for idx_cluster =  1:n_clusters
    
    [f, x] = ecdf(bouts_proc.durations_s(idx == idx_cluster));

    plot(x, f, 'Color', col.pca(idx_cluster, :), 'LineWidth', 2.5)

end
apply_generic(gca, 'xlim', [0 10.5], 'ylim', [0 1])
ylabel('ecdf')
xlabel('Time')
exporter(fh, paths, 'ecdf.pdf')
exporter(fh, paths, 'ecdf.png')