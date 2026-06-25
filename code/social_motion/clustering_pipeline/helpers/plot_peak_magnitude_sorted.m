function plot_peak_magnitude_sorted(al, clu, n_clusters, col, paths, varargin)
% PLOT_PEAK_MAGNITUDE_SORTED  4x3 figure: top row = break-density / mean traces
% (original, peak-aligned, magnitude-sorted), bottom = the matching heatmaps.
% Covers both the standard figure and the poster crop via options.
%
% Options (defaults = the standard figure):
%   position            figure Position
%   heatmap_clim        clim for the original & peak-aligned heatmaps
%   band_clusters       which clusters get a side colour band ([] -> 1:n_clusters)
%   yline_clusters      which boundaries get a dashed line  ([] -> 1:n_clusters)
%   show_cluster_means  overlay per-cluster mean traces on the peak panel
%   peak_alpha          mask NaNs on the peak-aligned heatmap (AlphaData)
%   export_name         '' to skip export, else base filename (no extension)

opt = inputParser;
addParameter(opt, 'position', [100 100 1800 1000]);
addParameter(opt, 'heatmap_clim', [0 3.2]);
addParameter(opt, 'band_clusters', []);
addParameter(opt, 'yline_clusters', []);
addParameter(opt, 'show_cluster_means', true);
addParameter(opt, 'peak_alpha', false);
addParameter(opt, 'export_name', '');
parse(opt, varargin{:});
o = opt.Results;

band_clusters  = o.band_clusters;  if isempty(band_clusters),  band_clusters  = 1:n_clusters; end
yline_clusters = o.yline_clusters; if isempty(yline_clusters), yline_clusters = 1:n_clusters; end

% unpack compute-stage outputs
boundaries             = clu.boundaries;
temp_cluster_ids       = clu.temp_cluster_ids;
break_times_valid      = al.break_times_valid;
break_times_peak_valid = al.break_times_peak_valid;
break_times_peak_mag   = al.break_times_peak_mag;
sorted_matrix_valid    = al.sorted_matrix_valid;
aligned_peak_valid     = al.aligned_peak_valid;
aligned_peak_mag       = al.aligned_peak_mag;
x_axis_peak            = al.x_axis_peak;
n_cols                 = al.n_cols;
n_rows_valid           = al.n_rows_valid;
n_rows_mag             = al.n_rows_mag;
n_mag_clusters         = al.n_mag_clusters;
idx_mag_cluster_sorted = al.idx_mag_cluster_sorted;
boundaries_sm          = al.boundaries_sm;

new_boundaries = boundaries - boundaries(2);
new_boundaries(1) = [];
new_boundaries(end) = [];

fh = figure('color','w','Position',o.position);
tiledlayout(4, 3, 'TileSpacing', 'compact', 'Padding', 'compact')
colsm = cbrewer2('Reds', 1);

bin_size = 10;

% ================= TOP ROW =================

% --- ORIGINAL
nexttile
ax_distr(1) = gca;
hold on
histogram(break_times_valid, 0:bin_size:90000, 'Normalization','pdf', ...
    'FaceColor','k','EdgeColor','none')

apply_generic(ax_distr(1), 'ylim', [0 .02], 'xlim', [0 630], 'font_size', 18, 'yticks', [0 0.01 0.02]);
ylabel('Break Density')

yyaxis right
avg = mean(sorted_matrix_valid,1,'omitnan');
plot(1:n_cols, avg, 'r-', 'LineWidth', 2, 'Color', colsm)
xline(1,'k--')
ylim([0 3])
set(gca,'XTick',[], 'FontSize', 18)
ylabel('Mean SM')
ax_distr(1).YAxis(2).Color = colsm;

% --- PEAK ALIGNED
nexttile
ax_distr(2) = gca;
hold on

edges = -601:bin_size:max(break_times_peak_valid);
[counts, ~] = histcounts(break_times_peak_valid, edges);

available_per_bin = sum(~isnan(aligned_peak_valid),1);
available_binned = zeros(size(counts));
for i = 1:length(counts)
    bin_mask = x_axis_peak >= edges(i) & x_axis_peak < edges(i+1);
    available_binned(i) = sum(available_per_bin(bin_mask));
end

counts_corrected = counts ./ available_binned;
counts_corrected = counts_corrected ./ (sum(counts_corrected) * bin_size);

histogram('BinEdges', edges, 'BinCounts', counts_corrected, ...
    'FaceColor', 'k', 'EdgeColor', 'none');

apply_generic(ax_distr(2), 'ylim', [0 .005], 'xlim', [-630 630]./2, 'font_size', 18, 'yticks', [0 0.01 0.02]);
ylabel('Break Density')

yyaxis right
if o.show_cluster_means
    for idx_valid_cluster = 1:n_clusters
        avg = mean(aligned_peak_valid(temp_cluster_ids == idx_valid_cluster,:), 1,'omitnan');
        plot(x_axis_peak, avg, 'r-', 'LineWidth', 1.5, 'Color', col.pca(idx_valid_cluster + 1 , :))
    end
end

avg = mean(aligned_peak_valid,1,'omitnan');
plot(x_axis_peak, avg, 'k-', 'LineWidth', 2, 'Color', colsm)
ylabel('Mean SM')

xline(1,'k--')
ylim([0 3])
set(gca,'XTick',[], 'FontSize', 18)
ax_distr(2).YAxis(2).Color = colsm;

% --- MAGNITUDE SORTED
nexttile
ax_distr(3) = gca;
hold on
colmag = flipud(cbrewer2('Set2', n_mag_clusters));

for k = 1:n_mag_clusters
    data = break_times_peak_mag(idx_mag_cluster_sorted == k);
    data = data(~isnan(data));
    [f, x] = ecdf(data);
    plot(x, f, 'Color', colmag(k,:), 'LineWidth', 2);
end

ylabel('CDF')
set(gca,'FontSize',18)
apply_generic(ax_distr(3), 'ylim', [0 1], 'xlim', [-630 630]./2, 'font_size', 18)
yyaxis right
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
imagesc(1:n_cols, 1:n_rows_valid, sorted_matrix_valid, o.heatmap_clim)
set(gca,'YDir','normal')
hold on
scatter(break_times_valid, 1:n_rows_valid, 3, '|', 'k')
apply_generic(gca, 'no_y', true, 'xlim', [0 630], 'xticks', 0:120:600)
xticklabels(0:2:10)
colormap(cbrewer2('Reds',[]))

for idx_cluster = band_clusters
    fill([ax_distr(4).XLim(1)-2,ax_distr(4).XLim(1)-2,ax_distr(4).XLim(1)-17,ax_distr(4).XLim(1)-17], [new_boundaries(idx_cluster)  new_boundaries(idx_cluster +1) new_boundaries(idx_cluster +1)   new_boundaries(idx_cluster)], ...
        col.pca(idx_cluster,:), 'EdgeColor','none', 'Clipping', 'off');
    text(mean([ax_distr(4).XLim(1)-2,ax_distr(4).XLim(1)-17]), mean([new_boundaries(idx_cluster); new_boundaries(idx_cluster + 1)]), num2str(idx_cluster + 1),...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
    hold on
end

xlabel('Time (s)');

hold on;
for idx_cluster = yline_clusters
    yline(new_boundaries(idx_cluster), 'k--', 'LineWidth', 1.1);
end

% --- PEAK ALIGNED
nexttile(5,[3 1])
ax_distr(5) = gca;
h = imagesc(x_axis_peak, 1:n_rows_valid, aligned_peak_valid, o.heatmap_clim);
if o.peak_alpha
    set(h, 'AlphaData', ~isnan(aligned_peak_valid));
end
set(gca,'YDir','normal')
hold on
scatter(break_times_peak_valid, 1:n_rows_valid, 3, '|', 'k')

apply_generic(gca, 'no_y', true, 'xlim', [-630 630]./2, 'xticks', -360:120:360)
xticklabels(-6:2:6)
colormap(cbrewer2('Reds',[]))

for idx_cluster = band_clusters
    fill([ax_distr(5).XLim(1)-2,ax_distr(5).XLim(1)-2,ax_distr(5).XLim(1)-17,ax_distr(5).XLim(1)-17], [new_boundaries(idx_cluster)  new_boundaries(idx_cluster+1) new_boundaries(idx_cluster+1)   new_boundaries(idx_cluster)], ...
        col.pca(idx_cluster,:), 'EdgeColor','none', 'Clipping', 'off');
    text(mean([ax_distr(5).XLim(1)-2,ax_distr(5).XLim(1)-17]), mean([new_boundaries(idx_cluster ); new_boundaries(idx_cluster +1)]), num2str(idx_cluster +1 ),...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
    hold on
end

xlabel('Time (s)');

hold on;
for idx_cluster = yline_clusters
    yline(new_boundaries(idx_cluster), 'k--', 'LineWidth', 1.1);
end

% --- MAG SORTED
nexttile(6,[3 1])
ax_distr(6) = gca;
imagesc(x_axis_peak, 1:n_rows_mag, aligned_peak_mag, [0 3.2])
set(gca,'YDir','normal')
hold on
scatter(break_times_peak_mag, 1:n_rows_mag, 3, '|', 'k')
apply_generic(gca, 'no_y', true, 'xlim', [-630 630]./2, 'xticks', -360:120:360)
xticklabels(-6:2:6)
colormap(cbrewer2('Reds',[]))

x_left = ax_distr(6).XLim(1);
for k = 1:n_mag_clusters
    y1 = boundaries_sm(k) + 1;
    y2 = boundaries_sm(k + 1);
    fill([x_left-2, x_left-2, x_left-17, x_left-17], [y1 y2 y2 y1], colmag(k,:), 'EdgeColor','none', 'Clipping','off');
    text(x_left-9.5, mean([y1 y2]), num2str(k), 'HorizontalAlignment','center', 'VerticalAlignment','middle');
end

xlabel('Time (s)');

% ================= LINK =================
linkaxes([ax_distr(1) ax_distr(4)], 'x')
linkaxes([ax_distr(2) ax_distr(5)], 'x')
linkaxes([ax_distr(3) ax_distr(6)], 'x')

set(gcf, 'GraphicsSmoothing', 'off');

if ~isempty(o.export_name)
    exporter(fh, paths, [o.export_name '.pdf'])
    exporter(fh, paths, [o.export_name '.png'])
end
end
