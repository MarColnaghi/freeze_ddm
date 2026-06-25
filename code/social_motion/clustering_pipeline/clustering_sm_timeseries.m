
clearvars

% Stage 1 — load loom-evoked freezes + pre-freeze social motion (contact-masked).
% See load_prefreeze_sm for the options (le_window, min_dur, nloom, ...).
data = load_prefreeze_sm();
bouts_proc         = data.bouts_proc;
sm_freeze_full     = data.sm_freeze_full;   % trials x time
is_below_threshold = data.is_below_threshold;
inizio             = data.inizio;

fine               = data.fine;

% sm_freeze_full is trials x timesteps (each row is a trial, each column a
% timestep in the duration)

% Figure-output paths (per clustering type)
clustering_type = 'max';
threshold_imm = 2; threshold_mob = 2; threshold_pc = 4; id_code = sprintf('imm%d_mob%d_pc%d', threshold_imm, threshold_mob, threshold_pc);
paths = path_generator('folder', fullfile('social_motion','clustering', clustering_type), 'bouts_id', id_code, 'imfirst', false);

%% 
% Visualize the Heatmap
fh = figure('color','w','Position',[100, 100, 600, 4000]);
tiledlayout(4, 1, 'Padding', 'compact', 'TileSpacing', 'tight')
nexttile(2, [3 1])
ax_distr(2) = gca;
[~, sort_order] = sort(bouts_proc.sm);

hold on
x_axis = inizio:fine;
y_axis = 1:size(sm_freeze_full(sort_order, :), 1);
h = imagesc(x_axis, y_axis, sm_freeze_full(sort_order, :), [0 2]);
% set(h, 'AlphaData', ~isnan(sm_freeze(sort_order, :)));

scatter(bouts_proc.durations(sort_order), 1:height(bouts_proc), 2, '|', 'k')
apply_generic(gca, 'xlim', [0 630], 'no_y', true, 'ylim', [-100 size(sm_freeze_full(sort_order, :), 1) + 100], 'xtick', 0:120:600)
colormap(cbrewer2('Reds',[]));

%cb = colorbar(ax, 'Location', 'southoutside', 'FontSize', 18, 'LineWidth', 2);
%cb.Label.String = 'Social Motion';
xticklabels({'Freeze\newlineOnset', '2', '4', '6', '8', '10'});

% Set the interpreter to TeX so it recognizes the \newline command
ax = gca;
ax.TickLabelInterpreter = 'tex';
xtickangle(0);

xlabel('Time (s)');

nexttile(1)
ax_distr(1) = gca;
hold on
apply_generic(gca, 'ylim', [0 1.2],'ytick', [0 1.2], 'xlim', [0 630], 'no_x', true, 'font_size', 18, 'xticks', 0:120:600);
ylabel({'Mean', 'Social Motion'})

linkaxes([ax_distr(:)], 'x')
axi = gca;
axi.YAxis.Color = [1 1 1];
axi.XAxis.Color = [1 1 1];

exporter(fh, paths, 'notordered_profiles.pdf')
exporter(fh, paths, 'notordered_profiles.png')

%%
n_clusters = 10;

% Stage 2 — cluster the pre-freeze SM (see cluster_sm)
clu = cluster_sm(sm_freeze_full, bouts_proc, clustering_type, n_clusters);
idx                 = clu.idx;
sort_order          = clu.sort_order;
valid_mask          = clu.valid_mask;
invalid_mask        = clu.invalid_mask;
temp_cluster_ids    = clu.temp_cluster_ids;
repr                = clu.repr;
all_unique_clusters = clu.all_unique_clusters;
boundaries          = clu.boundaries;
sorted_matrix       = clu.sorted_matrix;

%% 

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
h = imagesc(x_axis, y_axis, sm_freeze_full(sort_order, :), [0 2]);
% set(h, 'AlphaData', ~isnan(sm_freeze(sort_order, :)));

scatter(bouts_proc.durations(sort_order), 1:height(bouts_proc), 2, '|', 'k')
apply_generic(gca, 'xlim', [0 630], 'no_y', true, 'ylim', [-100 size(sm_freeze_full, 1) + 100], 'xtick', 0:120:600)
colormap(cbrewer2('Reds',[]));

%cb = colorbar(ax, 'Location', 'southoutside', 'FontSize', 18, 'LineWidth', 2);
%cb.Label.String = 'Social Motion';
xticklabels({'Freeze\newlineOnset', '2', '4', '6', '8', '10'});

% Set the interpreter to TeX so it recognizes the \newline command
ax = gca;
ax.TickLabelInterpreter = 'tex';
xtickangle(0);

xlabel('Time (s)');

%cb = colorbar(ax, 'Location', 'southoutside', 'FontSize', 18, 'LineWidth', 2);
%cb.Label.String = 'Social Motion';

for idx_cluster = 1:length(all_unique_clusters)
    fill([ax_distr(2).XLim(1)-2,ax_distr(2).XLim(1)-2,ax_distr(2).XLim(1)-17,ax_distr(2).XLim(1)-17], [boundaries(idx_cluster) boundaries(idx_cluster + 1) boundaries(idx_cluster + 1)   boundaries(idx_cluster)], ...
        col.pca(idx_cluster,:), 'EdgeColor','none', 'Clipping', 'off');
    text(mean([ax_distr(2).XLim(1)-2,ax_distr(2).XLim(1)-17]), mean([boundaries(idx_cluster); boundaries(idx_cluster + 1)]), num2str(idx_cluster),...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
    hold on

end

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
xlim([0 630]);

exporter(fh, paths, 'clusters_profiles.pdf')
exporter(fh, paths, 'clusters_profiles.png')

axes(ax_distr(2))
clim([90 91.2])
exporter(fh, paths, 'clusters_profiles_nocolor.pdf')
exporter(fh, paths, 'clusters_profiles_nocolor.png')

axes(ax_distr(2))
clim([0 2])
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


%% --- Peak/fit alignment (see align_peaks)
al = align_peaks(clu, sm_freeze_full, bouts_proc, inizio);
valid_sorted           = al.valid_sorted;
x_axis_peak            = al.x_axis_peak;
aligned_peak_valid     = al.aligned_peak_valid;
aligned_peak_mag       = al.aligned_peak_mag;
sorted_matrix_valid    = al.sorted_matrix_valid;
break_times_valid      = al.break_times_valid;
break_times_peak_valid = al.break_times_peak_valid;
break_times_peak_mag   = al.break_times_peak_mag;
n_rows_valid           = al.n_rows_valid;
n_rows_mag             = al.n_rows_mag;
n_cols                 = al.n_cols;
n_mag_clusters         = al.n_mag_clusters;
idx_mag_cluster_sorted = al.idx_mag_cluster_sorted;
boundaries_sm          = al.boundaries_sm;
%% --- Peak-vs-magnitude figure (standard)
plot_peak_magnitude_sorted(al, clu, n_clusters, col, paths, ...
    'export_name', 'peak_vs_magnitude_sorted');
%% --- Peak-vs-magnitude figure (poster crop)
plot_peak_magnitude_sorted(al, clu, n_clusters, col, paths, ...
    'position', [100 100 1200 700], 'heatmap_clim', [0 5], ...
    'band_clusters', 3, 'yline_clusters', 3:4, ...
    'show_cluster_means', false, 'peak_alpha', true, ...
    'export_name', 'peak_vs_magnitude_sorted_4poster');
%%
colors = cbrewer2('Set1', 3);
colors = colors([2, 1, 3], :);

% SM aligned to freeze offset (before/during/after on peak; duration split on offset)
sm_freeze = extract_sm_from_bouts(bouts_proc, 'type', 'onlyfreeze', 'output_type', 'mat', 'align', 'offset', 'delay_after', 60, 'norm_factor', 10);
plot_aligned_offset(sm_freeze, al, sort_order, bouts_proc, colors, paths, 'sm_aligned_2_offset', true);


%%


n_trials = 400;
selected_trials = sort(randi(size(aligned_peak_valid, 1), n_trials, 1));

fh = figure('color','w','Position',[100 100 400 1200]);
tiledlayout(1, 1, 'TileSpacing', 'tight', 'Padding', 'tight')
nexttile
hold on
for idx_trial = 1:n_trials
    n = selected_trials(idx_trial);
    plot(x_axis_peak, aligned_peak_valid(n, :) + idx_trial, 'k-')
end

% x positions
x_breaks = break_times_peak_valid(selected_trials);

% row indices
rows = selected_trials;

% convert break times to column indices
cols = x_breaks + 1;

% extract corresponding y values
% scatter
apply_generic(gca, 'no_y', true, 'xlim', [-600 600]./2, 'xticks', -360:120:360, 'no_x', true)
xticklabels(-6:2:6)

exporter(fh, paths, 'sm_traces_aligned_2_peak.pdf')
exporter(fh, paths, 'sm_traces_aligned_2_peak.png')

%%
fh = figure('Position', [100 100 1300 700], 'Color', 'w');



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

n_clusters_grid = numel(all_unique_clusters);
n_cols_grid = 4;
n_rows_grid = ceil(n_clusters_grid / n_cols_grid);

tiledlayout(n_rows_grid, n_cols_grid, 'Padding','loose','TileSpacing','compact');

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

for idx_cluster = 1:n_clusters_grid

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

% Cluster band boundaries (the peak-vs-magnitude cells used to leave this in the
% workspace; now that they are functions, compute it here).
new_boundaries = boundaries - boundaries(2);
new_boundaries(1) = [];
new_boundaries(end) = [];

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


%% Distance aligned to freeze offset (same panels, on min inter-fly distance)
distance = extract_sm_from_bouts(bouts_proc, 'type', 'onlyfreeze', 'output_type', 'mat', 'align', 'offset', 'delay_after', 60, 'cache', 'mindist_cache');
plot_aligned_offset(distance, al, sort_order, bouts_proc, colors, paths, '', false);
