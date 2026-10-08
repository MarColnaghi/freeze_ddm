id_code = 'imm3_mob3_pc4';
paths_out = path_generator('folder', 'spatial', 'bouts_id', id_code, 'imfirst', false);
if ~isfolder(paths_out.fig), mkdir(paths_out.fig); end

col = cmapper();
thresholds = define_thresholds;

bouts = importdata(fullfile(paths_out.dataset, 'bouts.mat'));
bouts = bouts_formatting(bouts, thresholds);

%% Map every fly onto the unit disk of its own arena
% The arena sits at a different place in the camera frame from one recording
% to the next (centres spread over ~10 mm), and its radius in mm differs too
% (~22-26 mm, tracking the position, so likely two rigs or calibrations).
% Pooling mm coordinates smears the wall, so each fly is re-centred on its
% arena and scaled by its radius, both fitted once by build_arena_cache from
% every position of all five flies in the raw files.

arena_file = fullfile(paths_out.cache_path, 'arena_cache.mat');
if ~isfile(arena_file)
    build_arena_cache();
end
arenas = importdata(arena_file);

[G, fly_ids] = findgroups(bouts.fly);
[~, idx_arena] = ismember(fly_ids, arenas.fly);
assert(all(idx_arena > 0), 'Some flies in bouts have no fitted arena.');
centre = [arenas.centre_x(idx_arena), arenas.centre_y(idx_arena)];
radius = arenas.radius(idx_arena);

% The fitted circle runs through the middle of the outermost positions, so
% flies hugging the wall land just past r = 1. Those are pulled back onto the
% wall, so everything lies in [-1, 1]
to_unit = @(xy) xy ./ max(1, hypot(xy(:, 1), xy(:, 2)));

xy_n = ([bouts.x_onset, bouts.y_onset] - centre(G, :)) ./ radius(G);
r_n = hypot(xy_n(:, 1), xy_n(:, 2));
fprintf('%.1f%% of bout onsets past the fitted wall, 99.9th percentile at r/R = %.3f\n', ...
    100 * mean(r_n > 1), prctile(r_n, 99.9));
xy_n = to_unit(xy_n);
bouts.x_c = xy_n(:, 1);
bouts.y_c = xy_n(:, 2);

arena_r = 1;

%% Loom-evoked freezes
min_dur = 30; % frames, 0 keeps every le immobility bout
bouts_le = data_parser_new(bouts, 'type', 'immobility', 'period', 'loom', 'window', 'le', 'min_dur', min_dur);

%% 2D density of onset positions, per loom speed
% Densities are in units of a uniform spread over the arena, so 1 means as
% many onsets per unit area as if freezes landed anywhere with equal probability
bw = 0.04;                               % kernel bandwidth (r/R, ~1 mm)
grid_xy = linspace(-arena_r, arena_r, 151);
dx = grid_xy(2) - grid_xy(1);
[xq, yq] = meshgrid(grid_xy);
outside = hypot(xq, yq) > arena_r;
arena_area = pi * arena_r^2;

% Radial profile: onsets per annulus, divided by the annulus area
r_edges = linspace(0, arena_r, 26);
r_centres = (r_edges(1:end-1) + r_edges(2:end)) / 2;

speeds = [25 50];
speed_names = {'Slow loom', 'Fast loom'};
speed_cols = col.vars.sloom([3 5], :);
dens = cell(1, numel(speeds));
dens_r = zeros(numel(speeds), numel(r_centres));
n_bouts = zeros(1, numel(speeds));
n_flies = zeros(1, numel(speeds));

for k = 1:numel(speeds)
    b = bouts_le(bouts_le.sloom == speeds(k), :);
    n_bouts(k) = height(b);
    n_flies(k) = numel(unique(b.fly));

    f = ksdensity([b.x_c, b.y_c], [xq(:), yq(:)], 'Bandwidth', [bw bw]);
    f = reshape(f, size(xq));
    f(outside) = NaN;
    % Renormalise inside the arena, the kernel leaks mass past the wall
    dens{k} = f / (sum(f(:), 'omitnan') * dx^2) * arena_area;

    r = hypot(b.x_c, b.y_c);
    dens_r(k, :) = histcounts(r, r_edges) / n_bouts(k) ./ (pi * diff(r_edges.^2)) * arena_area;

    fprintf('%s (%d): n=%d freezes, %d flies, median r/R=%.2f, within 0.2 R of wall=%.0f%%\n', ...
        speed_names{k}, speeds(k), n_bouts(k), n_flies(k), median(r), 100*mean(r > arena_r - 0.2));
end

dens_max = max(cellfun(@(f) max(f(:)), dens));

%% At-risk occupancy during the le windows
% Flies sit at the wall most of the time, so the onset density above mostly
% tracks where they are. The baseline here is where a freeze could start:
% every frame inside a loom's le window (same per-speed window as bouts.le) on
% which the fly is mobile. Bouts tile every frame, so the bout types repeated
% over their durations give a per-frame immobility mask.
xy_cache = importdata(fullfile(paths_out.cache_path, 'flyxy_cache.mat'));
loom_cache = importdata(fullfile(paths_out.cache_path, 'loom_cache.mat'));

occ_xy = cell(numel(fly_ids), 1);
fly_sloom = zeros(numel(fly_ids), 1);
for k = 1:numel(fly_ids)
    bf = sortrows(bouts(G == k, :), 'onsets');
    fly_sloom(k) = bf.sloom(1);
    if fly_sloom(k) == 25
        w = thresholds.le_window_sl;
    else
        w = thresholds.le_window_fl;
    end

    immobile = repelem(bf.type, bf.durations);
    looming = loom_cache(fly_ids(k));
    loom_starts = find(diff([0; looming(:)]) == 1);
    frames = loom_starts + (w(1):w(2));
    frames = frames(frames <= numel(immobile));
    frames = frames(~immobile(frames));

    xy = xy_cache(fly_ids(k));
    occ_xy{k} = to_unit((xy(frames, :) - centre(k, :)) / radius(k));
end
clear xy_cache

%% Freeze rate relative to occupancy
% Onsets per at-risk frame, as a rate map: onset counts and occupancy are
% binned on the same grid and smoothed with the same Gaussian before dividing.
% Each speed is scaled by its own arena-wide rate, so 1 means no spatial
% preference and the two speeds compare on spatial modulation alone
bw_rate = 0.08;                           % smoothing sigma (r/R, ~2 mm)
min_occ = 0.1;                            % mask bins below this x the mean smoothed occupancy
n_boot = 1000;
edges_map = [grid_xy - dx/2, grid_xy(end) + dx/2];
% Equal-area annuli, so the sparse centre gets as much arena as the wall
r_edges_rate = arena_r * sqrt((0:10) / 10);
r_centres_rate = (r_edges_rate(1:end-1) + r_edges_rate(2:end)) / 2;

% Per-fly counts per annulus, for the fly bootstrap
n_r = numel(r_centres_rate);
onsets_r = zeros(numel(fly_ids), n_r);
occ_r = zeros(numel(fly_ids), n_r);
[~, fly_idx_le] = ismember(bouts_le.fly, fly_ids);
for k = 1:numel(fly_ids)
    occ_r(k, :) = histcounts(hypot(occ_xy{k}(:, 1), occ_xy{k}(:, 2)), r_edges_rate);
    b = bouts_le(fly_idx_le == k, :);
    onsets_r(k, :) = histcounts(hypot(b.x_c, b.y_c), r_edges_rate);
end

rate_map = cell(1, numel(speeds));
rate_r = zeros(numel(speeds), n_r);
rate_ci = zeros(2, n_r, numel(speeds));
rng(1)

for k = 1:numel(speeds)
    b = bouts_le(bouts_le.sloom == speeds(k), :);
    occ = vertcat(occ_xy{fly_sloom == speeds(k)});
    mean_rate = height(b) / size(occ, 1);

    on_s = imgaussfilt(histcounts2(b.y_c, b.x_c, edges_map, edges_map), bw_rate / dx);
    occ_s = imgaussfilt(histcounts2(occ(:, 2), occ(:, 1), edges_map, edges_map), bw_rate / dx);
    f = on_s ./ occ_s / mean_rate;
    f(outside | occ_s < min_occ * mean(occ_s(~outside))) = NaN;
    rate_map{k} = f;

    O = onsets_r(fly_sloom == speeds(k), :);
    C = occ_r(fly_sloom == speeds(k), :);
    rate_r(k, :) = sum(O) ./ sum(C) / (sum(O(:)) / sum(C(:)));

    boot = zeros(n_boot, n_r);
    for i = 1:n_boot
        s = randi(size(O, 1), size(O, 1), 1);
        o = sum(O(s, :));
        c = sum(C(s, :));
        boot(i, :) = o ./ c / (sum(o) / sum(c));
    end
    rate_ci(:, :, k) = prctile(boot, [2.5 97.5]);

    fprintf('%s (%d): %.2f le freezes per at-risk second, %.1f h at risk\n', ...
        speed_names{k}, speeds(k), 60 * mean_rate, size(occ, 1) / 60 / 3600);
end

% Diverging scale on log2(rate), symmetric about 1x, saturating the top 1%
rate_lim = ceil(max(cellfun(@(f) prctile(abs(log2(f(:))), 99), rate_map)) * 2) / 2;

%% Plot
fh = figure('color', 'w', 'Position', [100, 100, 1500, 1040]);
tl = tiledlayout(2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

theta = linspace(0, 2*pi, 400);
cmap_dens = cbrewer2('Blues', 256);
map_lim = [-1 1] * 1.05 * arena_r;

for k = 1:numel(speeds)
    ax = nexttile;
    hold on

    f = dens{k};
    imagesc(grid_xy, grid_xy, f, 'AlphaData', ~isnan(f))
    colormap(ax, cmap_dens)
    clim(ax, [0 dens_max])
    plot(arena_r * cos(theta), arena_r * sin(theta), 'Color', [.3 .3 .3], 'LineWidth', 1.5)

    axis(ax, 'equal')
    xlim(ax, map_lim)
    ylim(ax, map_lim)
    apply_generic(ax, 'font_size', 20, 'line_width', 2, 'xticks', [-1 0 1], 'yticks', [-1 0 1])
    title(sprintf('%s (%d)\nn = %d freezes, %d flies', speed_names{k}, speeds(k), n_bouts(k), n_flies(k)), ...
        'FontWeight', 'normal', 'FontSize', 20)

    if k == numel(speeds)
        % Shared colour scale for both maps
        cb = colorbar(ax, 'eastoutside');
        cb.FontSize = 16;
        cb.LineWidth = 1.5;
        cb.Label.String = 'Density (\times uniform)';
    end
end

ax = nexttile;
hold on
yline(1, ':', 'Color', [.5 .5 .5], 'LineWidth', 1.5)
lh = gobjects(1, numel(speeds));
for k = 1:numel(speeds)
    lh(k) = plot(r_centres, dens_r(k, :), 'Color', speed_cols(k, :), 'LineWidth', 2.5, ...
        'DisplayName', sprintf('%s (%d)', speed_names{k}, speeds(k)));
end
axis(ax, 'square')
xlim(ax, [0 arena_r])
ylim(ax, [0 1.1 * max(dens_r(:))])
apply_generic(ax, 'font_size', 20, 'line_width', 2)
xlabel('Distance from centre (r / R)')
ylabel('Density (\times uniform)')
legend(ax, lh, 'Location', 'northwest', 'Box', 'off', 'FontSize', 16)
title('Radial profile', 'FontWeight', 'normal', 'FontSize', 20)

% Second row: onset rate relative to at-risk occupancy
cmap_rate = flipud(cbrewer2('RdBu', 256));
rate_ticks = 2 .^ (-floor(rate_lim):floor(rate_lim));

for k = 1:numel(speeds)
    ax = nexttile;
    hold on

    % Grey underlay shows through the masked low-occupancy bins
    fill(arena_r * cos(theta), arena_r * sin(theta), [.85 .85 .85], 'EdgeColor', 'none')
    f = log2(rate_map{k});
    imagesc(grid_xy, grid_xy, f, 'AlphaData', ~isnan(f))
    colormap(ax, cmap_rate)
    clim(ax, [-1 1] * rate_lim)
    plot(arena_r * cos(theta), arena_r * sin(theta), 'Color', [.3 .3 .3], 'LineWidth', 1.5)

    axis(ax, 'equal')
    xlim(ax, map_lim)
    ylim(ax, map_lim)
    apply_generic(ax, 'font_size', 20, 'line_width', 2, 'xticks', [-1 0 1], 'yticks', [-1 0 1])
    title(sprintf('%s (%d)\nrelative to occupancy', speed_names{k}, speeds(k)), ...
        'FontWeight', 'normal', 'FontSize', 20)

    if k == numel(speeds)
        cb = colorbar(ax, 'eastoutside');
        cb.FontSize = 16;
        cb.LineWidth = 1.5;
        cb.Ticks = log2(rate_ticks);
        cb.TickLabels = compose('%g', rate_ticks);
        cb.Label.String = 'Freeze rate (\times mean)';
    end
end

ax = nexttile;
hold on
yline(1, ':', 'Color', [.5 .5 .5], 'LineWidth', 1.5)
lh = gobjects(1, numel(speeds));
for k = 1:numel(speeds)
    fill([r_centres_rate, fliplr(r_centres_rate)], [rate_ci(1, :, k), fliplr(rate_ci(2, :, k))], ...
        speed_cols(k, :), 'EdgeColor', 'none', 'FaceAlpha', 0.25)
    lh(k) = plot(r_centres_rate, rate_r(k, :), 'Color', speed_cols(k, :), 'LineWidth', 2.5, ...
        'DisplayName', sprintf('%s (%d)', speed_names{k}, speeds(k)));
end
axis(ax, 'square')
xlim(ax, [0 arena_r])
ylim(ax, [0 1.1 * max(rate_ci(:))])
apply_generic(ax, 'font_size', 20, 'line_width', 2)
xlabel('Distance from centre (r / R)')
ylabel('Freeze rate (\times mean)')
legend(ax, lh, 'Location', 'northwest', 'Box', 'off', 'FontSize', 16)
title(sprintf('Radial, 95%% CI over flies'), 'FontWeight', 'normal', 'FontSize', 20)

exporter(fh, paths_out, 'le_onset_xy_density.pdf')
