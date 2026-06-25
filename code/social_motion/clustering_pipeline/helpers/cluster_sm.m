function clu = cluster_sm(sm_freeze_full, bouts_proc, clustering_type, n_clusters)
% CLUSTER_SM  Cluster pre-freeze social-motion timeseries and order trials.
%
%   clu = cluster_sm(sm_freeze_full, bouts_proc, clustering_type, n_clusters)
%
% clustering_type: 'max' | 'com' | 'pca' | 'first_cross' (only 'max'/'com' are
% exercised end-to-end downstream).
%
% Fields of CLU (the outputs consumed by the plotting/alignment stages):
%   idx                 cluster id per trial
%   sort_order          row order (cluster, then metric)
%   valid_mask          trials with a usable peak ('max'/'com' only)
%   invalid_mask        trials that never reached threshold ('max'/'com' only)
%   temp_cluster_ids    within-valid cluster ids (length = #valid trials)
%   repr                cluster-mean profiles (n_clusters_present x time)
%   all_unique_clusters present cluster ids
%   boundaries          row boundaries between clusters (for plotting)
%   sorted_matrix       centered_sm_freeze reordered by sort_order

% Optional outputs only the 'max'/'com' paths define
valid_mask = []; invalid_mask = []; temp_cluster_ids = [];

switch clustering_type
    case 'pca'
        % 1. Run PCA
        centered_sm_freeze = zscore(sm_freeze_full, [], 2);
        % centered_sm_freeze = sm_freeze_full;

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
            threshold_sm = .5;
            [max_vals, metric] = max(centered_sm_freeze, [], 2);
            metric = double(metric);

            failed_to_reach = max_vals < threshold_sm;
            trunc_point = 30;
            early_peaks = metric < trunc_point;

        else
            % Metric: Center of Mass
            pow = 9;
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
        threshold = 0.75;
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
bouts_sorted = bouts_proc(sort_order, :); %#ok<NASGU>

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

clu = struct('idx', idx, 'sort_order', sort_order, ...
    'valid_mask', valid_mask, 'invalid_mask', invalid_mask, ...
    'temp_cluster_ids', temp_cluster_ids, 'repr', repr, ...
    'all_unique_clusters', all_unique_clusters, 'boundaries', boundaries, ...
    'sorted_matrix', sorted_matrix);
end
