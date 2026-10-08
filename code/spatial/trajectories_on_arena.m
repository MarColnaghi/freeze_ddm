paths_out = path_generator('folder', 'spatial/trajectories');
if ~isfolder(paths_out.fig), mkdir(paths_out.fig); end

% The position cache holds every fly and takes ~6 s to load, so it is loaded
% once and kept in the workspace across reruns
if ~exist('xys', 'var')
    xys = importdata(fullfile(paths_out.cache_path, 'flyxy_cache.mat'));
end
arenas = importdata(fullfile(paths_out.cache_path, 'arena_cache.mat'));

flies = arenas.fly';   % every fly, or a subset, e.g. [69 976]

% The fitted circle runs through the middle of the outermost positions, so
% flies hugging the wall land just past r = 1. Those are pulled back onto the
% wall, so everything lies in [-1, 1]
to_unit = @(xy) xy ./ max(1, hypot(xy(:, 1), xy(:, 2)));

arena_r = 1;
theta = linspace(0, 2*pi, 400);

% One hidden figure for the whole loop, only the trajectory changes between
% flies. A new figure per fly is ~3x slower
fh = figure('color', 'w', 'Position', [100, 100, 400, 400], 'Visible', 'off');
tl = tiledlayout(1, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
ax = nexttile;
hold on
plot(arena_r * cos(theta), arena_r * sin(theta), 'Color', [.3 .3 .3], 'LineWidth', 2)
traj = plot(nan, nan, 'LineWidth', 1.25, 'Color', '#DA4E76');
axis(ax, 'equal')
apply_generic(ax, 'font_size', 20, 'line_width', 2, 'xticks', [-1 0 1], 'yticks', [-1 0 1], ...
    'no_x', true, 'no_y', true, 'xlims', [-1.05 1.05], 'ylims', [-1.05 1.05])

for idx_fly = flies
    is_fly = arenas.fly == idx_fly;
    centre = [arenas.centre_x(is_fly), arenas.centre_y(is_fly)];
    radius = arenas.radius(is_fly) * 1.025;

    xy_norm = to_unit((xys(idx_fly) - centre) ./ radius);
    xy_norm = simplify_path(xy_norm, 0.005);

    set(traj, 'XData', xy_norm(:, 1), 'YData', xy_norm(:, 2))
    exporter(fh, paths_out, sprintf('fly_%03d.pdf', idx_fly))
end
close(fh)


function xy = simplify_path(xy, tol)
% Drops every point closer than tol to the last kept one, so no dropped point
% is further than tol from the drawn path. Frozen stretches collapse to a
% single vertex, which roughly halves the vertices and the PDF size

keep = false(size(xy, 1), 1);
keep(1) = true;
last = xy(1, :);
for j = 2:size(xy, 1)
    if hypot(xy(j, 1) - last(1), xy(j, 2) - last(2)) > tol
        keep(j) = true;
        last = xy(j, :);
    end
end
keep(end) = true;
xy = xy(keep, :);

end
