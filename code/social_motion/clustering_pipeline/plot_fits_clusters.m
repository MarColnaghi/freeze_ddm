model_list = {'dddim2'};
run_list   = {'run03_260407'};
version = {''};

paths_fit  = path_generator('folder', 'fitting_freezes/le_new');

paths_analysis = path_generator('folder', fullfile('social_motion', model_list{1}, strcat(run_list{1}, version{1})));
paths_analysis.fig = fullfile(paths_analysis.fig);
col = cmapper('', n_clusters);

model_results = importdata(fullfile(paths_fit.results, model_list{1}, strcat(run_list{1}, version{1}), 'model_results.mat'));
est_params.mean = model_results.estimates_mean;
est_params.std = model_results.estimates_std;  
extra = importdata(fullfile(paths_fit.results, model_list{1}, strcat(run_list{1}, version{1}), 'extra.mat'));


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
extra.soc_mot_array = extra.soc_mot_array(sort_order, :);

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

data = alternative_test_collect_data(model_list, run_list);

% Slice the tables to reduce communication overhead
out_all_tv = data(1).out_all_tv(sort_order, :);
out_all_st = data(1).out_all_st(sort_order, :);

% --- 3. Main Loop: Iterate over Bouts ---
n_bouts     = height(freezes_ref);
censoring_val = data(1).results.points.censoring;
res_points = data(1).results.points; % Slice to avoid passing large 'data' struct

fh = figure('color','w','Position',[100, 100, 1000, 550]);
tiledlayout(3, 4, 'Padding', 'compact', 'TileSpacing', 'compact')

fps = 60;
bin_size =  6 / fps;
t_vec = 0:1/fps:10.5;

for idx_cluster = 1:n_clusters
    
    bouts_cluster = bouts_sorted(cluster_id == idx_cluster, :);
    out_all_tv_cluster = out_all_tv(cluster_id == idx_cluster, :);
    out_all_st_cluster = out_all_st(cluster_id == idx_cluster, :);

    ax = nexttile;
    hold(ax, 'on')
    plot(t_vec, repr(idx_cluster, :), 'Color', col.pca(idx_cluster,:), 'LineWidth', 2)

    % Histogram
    histogram(ax, bouts_cluster.durations_s, ...
        min(bouts_proc.durations_s):bin_size:900, ...
        'Normalization', 'pdf', ...
        'FaceColor', col.pca(idx_cluster,:), ...
        'EdgeColor', 'k', ...
        'DisplayStyle', 'bar',...
        'LineWidth', 0.5);

    apply_generic(ax, 'xlim', [0 10], 'ylim', [-0.2 2.2])


    pdf_tv = zeros(631,1); pdf_tv_survival = 0;
    pdf_st = zeros(631,1); pdf_st_survival = 0;

    for idx_selected_bout = 1:height(bouts_cluster)
        
        fprintf('cluster: %d, bout: %d \n   ', idx_cluster, idx_selected_bout);

        % Current freeze row
        freeze_row = bouts_cluster(idx_selected_bout, :);
        dur_s      = freeze_row.durations_s;
        is_censored = dur_s > censoring_val;
        dur_s(is_censored) = censoring_val + 1/60;

        % Update external covariate for current bout
        ec.soc_mot_array = extra.soc_mot_array(idx_selected_bout, :);

        % --- Model 1 (Integration) ---
        cur_out_tv = out_all_tv_cluster(idx_selected_bout, :);
        [pdf_ddm_tv] = compute_pdf_tv_ddm(cur_out_tv, data(1).results.points);
        pdf_tv = pdf_tv + pdf_ddm_tv.ddm;
        pdf_tv_survival = pdf_tv_survival + pdf_ddm_tv.survival;

        cur_out_st = out_all_st_cluster(idx_selected_bout, :);
        [pdf_ddm_st] = compute_pdf_tv_ddm(cur_out_st, data(1).results.points);
        pdf_st = pdf_st + pdf_ddm_st.ddm;
        pdf_st_survival = pdf_st_survival + pdf_ddm_st.survival;

    end

    pdf_tv = pdf_tv ./ height(bouts_cluster);
    pdf_st = pdf_st ./ height(bouts_cluster);
    pdf_tv_survival = pdf_tv_survival ./ height(bouts_cluster);
    pdf_st_survival = pdf_st_survival ./ height(bouts_cluster);

    plot(pdf_ddm_tv.t, pdf_tv, 'k--', 'LineWidth', 2)
    plot(pdf_ddm_tv.t, pdf_st, 'r--', 'LineWidth', 2)

    %[~, fd, f] = nll_fly_ddm_newer(est_params.mean, bouts_sorted(cluster_id == idx_cluster))
 
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
    drawnow
end
