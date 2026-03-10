function [fh, ax] = fd_conditions_for_cosyne(varargin)
% --- Input Parsing ---
opt = inputParser;
addParameter(opt, 'extra', []);
addParameter(opt, 'freezes', []);
addParameter(opt, 'results', []);
addParameter(opt, 'bin_size', 3);
addParameter(opt, 'conditions', false);
addParameter(opt, 'type', 'continuous');
addParameter(opt, 'no_y', true);
addParameter(opt, 'color', 'col');
addParameter(opt, 'censored_inset', true);
addParameter(opt, 'gt', false);
addParameter(opt, 'export', false);
addParameter(opt, 'paths', []);
addParameter(opt, 'vis', 'off');
parse(opt, varargin{:});

results = opt.Results.results;
freezes = opt.Results.freezes;
censored_inset = opt.Results.censored_inset;
bin_size = opt.Results.bin_size;
color = opt.Results.color;
no_y = opt.Results.no_y;
hand_vis = opt.Results.vis;
export = opt.Results.export;
paths = opt.Results.paths;
conditions = opt.Results.conditions;

col = cmapper();
dx = 1/60; xxi = 0:dx:1200; h = .05; fs = 60;
bin_size_in_seconds = bin_size/fs;

if isempty(freezes)
    freezes = importdata(fullfile(results.bouts_path, 'surrogate.mat'));
end
if ~isempty(results.points.censoring)
    freezes.durations_s(freezes.durations_s > results.points.censoring) = results.points.censoring + 1;
end
freezes = freezes(freezes.durations_s >= results.points.truncation, :);

% --- Figure Setup ---
fh = figure('Position', [100 100 1200 900], 'Color', 'w');
t = tiledlayout(3, 4, 'TileSpacing', 'tight', 'Padding', 'compact');

visual_col_order = [1, 3, 2, 4];

i = 0;
for idx_sm = 1:3
    for j = 1:4
        i = i + 1;
        if j <= 2, idx_ls = 1; else, idx_ls = 2; end
        if mod(j,2) ~= 0, idx_fs = 1; else, idx_fs = 2; end

        target_col = visual_col_order(j);
        tile_idx = (idx_sm - 1) * 4 + target_col;

        ax(i) = nexttile(t, tile_idx);
        hold on;

        % --- Plotting Logic ---
        [freezes_quant, ~] = quantilizer_v2(freezes, 'indexed_quantile', ...
            struct('sm', idx_sm, 'ls', idx_ls, 'fs', idx_fs));

        censored = freezes_quant.durations_s > 10.5;
        [ks_noc] = ksdensity(freezes_quant.durations_s, xxi, 'BoundaryCorrection', 'reflection', 'Bandwidth', h, 'Support', [min(freezes.durations_s) - 0.001, max(freezes.durations_s) + 0.001 ]);
        [ks] = ksdensity(freezes_quant.durations_s, xxi, 'BoundaryCorrection', 'reflection', 'Bandwidth', h, 'Support', [min(freezes.durations_s) - 0.001, max(freezes.durations_s) + 0.001 ], 'Censoring', censored);
        censored_density = abs(trapz(xxi, ks) - trapz(xxi, ks_noc));

        histogram(freezes_quant.durations_s, min(freezes.durations_s) - 1/120:bin_size_in_seconds:900, 'Normalization', 'pdf', 'EdgeColor', 'none', 'FaceColor', col.empirical.(color)(2*(idx_sm - 1) + idx_fs, :), 'FaceAlpha', 0.25);
        plot(xxi, ks, 'Color', col.empirical.(color)(2*(idx_sm - 1) + idx_fs, :), 'LineWidth', 1.5);

        censored_x = results.points.censoring;
        scatter(censored_x, censored_density, 32, 'filled', 'MarkerFaceColor', col.empirical.(color)(2*(idx_sm - 1) + idx_fs, :), 'MarkerEdgeColor', 'none');

        apply_generic(ax(i), 'xlim', [-0.1 10.1], 'ylim', [-0.05 1.75], 'no_y', no_y, 'font_size', 20, 'line_width', 2.2, 'yticks', [0 1], 'xticks', [results.points.truncation 10], 'tick_length', 0.025);

        xlabel('Freeze Duration (s)')

        if i > 1
            xticklabels({});
            xlabel('')
        end
        % --- Inset Plot ---
        if censored_inset
            drawnow;
            cur_pos = ax(i).Position;
            inset_pos = [cur_pos(1) + 0.7*cur_pos(3), cur_pos(2) + 0.6*cur_pos(4), 0.02, 0.08];
            ax_inset(i) = axes('Position', inset_pos);
            hold(ax_inset(i), 'on'); set(ax_inset(i), 'box', 'on');
            scatter(ax_inset(i), censored_x, censored_density, 32, 'filled', 'MarkerFaceColor', col.empirical.(color)(2*(idx_sm - 1) + idx_fs, :), 'MarkerEdgeColor', 'none');
            apply_generic(ax_inset(i), 'ylim', [0 0.4], 'xlim', [censored_x - 0.025 censored_x + 0.025], 'yticks', [0 0.4], 'xticks', censored_x, 'tick_length', 0.05, 'line_width', 2, 'no_x', true, 'font_size', 12);
            if target_col > 1, yticklabels(ax_inset(i), {}); end
            axes(ax(i));
        end
    end
end

% --- ALIGNED COLOR PATCHES & TEXT (ASPECT RATIO FIXED) ---
drawnow;
fig_pos = fh.Position;
% For normalized units, we need the ratio of Width to Height in pixels
ratio = fig_pos(4) / fig_pos(3);

h_thick = 0.015; % Horizontal thickness
v_thick = h_thick * ratio; % Vertical thickness (scaled for pixels)

% --- 1. Top Spanning Boxes (Loom Type Indicators) ---
% Group 1: Columns 1 and 2
pos1 = []; pos2 = [];
% Group 2: Columns 3 and 4
pos3 = []; pos4 = [];

for k = 1:length(ax)
    if ax(k).Layout.Tile == 1, pos1 = ax(k).Position; end
    if ax(k).Layout.Tile == 2, pos2 = ax(k).Position; end
    if ax(k).Layout.Tile == 3, pos3 = ax(k).Position; end
    if ax(k).Layout.Tile == 4, pos4 = ax(k).Position; end
end

% Box for Col 1-2
if ~isempty(pos1) && ~isempty(pos2)
    span_width = (pos2(1) + pos2(3)) - pos1(1);
    annotation('rectangle', [pos1(1), pos1(2) + pos1(4) - 0, span_width, h_thick], ...
        'FaceColor', col.vars.fs(2,:), 'EdgeColor', 'none'); % Use specific color for Loom A
end

% Box for Col 3-4
if ~isempty(pos3) && ~isempty(pos4)
    span_width = (pos4(1) + pos4(3)) - pos3(1);
    annotation('rectangle', [pos3(1), pos3(2) + pos3(4) - 0, span_width, h_thick], ...
        'FaceColor', col.vars.fs(4,:), 'EdgeColor', 'none'); % Use specific color for Loom B
end

% --- 2. Standard Column Logic (Top Individual Patches) ---
col_labels = {'Slow Loom', 'Slow Loom', 'Fast Loom', 'Fast Loom'};
col_colors = {col.vars.ls(2,:), col.vars.ls(4,:), col.vars.ls(2,:), col.vars.ls(4,:)};

for c = 1:4
    target_handle = [];
    for k = 1:length(ax)
        if ax(k).Layout.Tile == c, target_handle = ax(k); break; end
    end
    if ~isempty(target_handle)
        pos = target_handle.Position;
        annotation('rectangle', [pos(1), pos(2) + pos(4) - 0.02, pos(3), h_thick], ...
            'FaceColor', col_colors{c}, 'EdgeColor', 'none');
    end
end

% --- 3. Row Logic (Left Individual Patches) ---
row_labels = {'Low Soc. Motion', 'Medium Soc. Motion', 'High Soc. Motion'};
row_colors = {col.vars.sm(1,:), col.vars.sm(2,:), col.vars.sm(3,:)};

for r = 1:3
    row_first_tile = (r-1)*4 + 1;
    target_handle = [];
    for k = 1:length(ax)
        if ax(k).Layout.Tile == row_first_tile, target_handle = ax(k); break; end
    end
    if ~isempty(target_handle)
        pos = target_handle.Position;
        % v_thick used here to match pixel-thickness of h_thick
        annotation('rectangle', [pos(1) - v_thick - 0.015, pos(2), v_thick, pos(4)], ...
            'FaceColor', row_colors{r}, 'EdgeColor', 'none');
    end
end
linkaxes(ax, 'xy');

if export
    paths.fig = results.fig_path;
    filename = 'fits.pdf'; if conditions, filename = 'fits_xcondition.pdf'; end
    exporter(fh, paths, filename);
end
end