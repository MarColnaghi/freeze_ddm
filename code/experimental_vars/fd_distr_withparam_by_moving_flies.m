function fh = fd_distr_withparam_by_moving_flies(varargin)
% FD_DISTR_WITHPARAM_BY_MOVING_FLIES  Freeze duration distribution against a control
% variable, one panel per number of moving flies.
%
%   fh = fd_distr_withparam_by_moving_flies('bouts', freezes, 'params', {'sm_1s'})
%
% The same figure as fd_distr_withparam_new, which has one panel per loom speed, with
% the panels split by the number of moving flies around the focal fly instead. Each
% panel shows the cumulative distribution of freeze duration at every level of the
% control variable (Kaplan-Meier, with bouts.is_censored as the censoring), and above
% it the distribution of the variable in that group, coloured by level. The colours
% are those of cmapper, and the circle of each panel is the colour of its number of
% moving flies.
%
% bouts is the table of freezes of loom_evoked_freezes, with durations_s, is_censored,
% ls, moving_flies and the control variables.
%
% Options
%   params           {'sm_1s', 'avg_fs_1s', 'avg_fs_loom', 'nloom'}   one figure per variable
%   type             'cumulative' (Kaplan-Meier), or 'kde', the Kaplan-Meier mass of
%                              every 1/3 s bin, or 'whisk' as in fd_distr_withparam_new.
%                              'whisk' does not take the censoring into account
%   n_quantiles      3         number of levels of a continuous variable, 3 for terciles.
%                              'nloom' and 'ls' keep their fixed levels
%   quantiles        'group'   the quantiles of a continuous variable are taken within
%                              each number of moving flies, so every panel compares the
%                              same number of freezes in each colour and the variable is
%                              compared within the same experimental group. The colours
%                              then mean different values in different panels.
%                    'pooled'  the quantiles of all the panels together, which is how
%                              fd_distr_withparam_new takes them
%   sloom_to_plot    'all'     'slow' or 'fast' keeps one loom speed, and the quantiles
%                              are those of the freezes that are plotted
%   moving_to_plot   'all'     or a vector of the numbers of moving flies, e.g. [0 4]
%   period           []        'loom' or 'bsl', taken from bouts.period if empty
%   export           false     writes fd_<variable>_by_moving_flies_<period>.pdf and .png
%                              to paths.fig
%   paths            []        struct with the figure folder in .fig, as exporter takes
%   check_quantiles  false     prints the edges and the number of freezes of every level
%
% 'nloom' and 'ls' have fixed levels (loom 2-5, 6-10, ... and slow, fast), the same in
% every panel. 'moving_flies' as the variable is not split, since the panels would be
% its levels: it gives one panel with a curve for every number of moving flies.
%
% A level with no freeze that ended before it was censored is left out of the panel.

%% ===================== Parse inputs =====================
opt = inputParser;
addParameter(opt, 'bouts', []);
addParameter(opt, 'export', false);
addParameter(opt, 'period', []);
addParameter(opt, 'paths', []);
addParameter(opt, 'params', {'sm_1s', 'avg_fs_1s', 'avg_fs_loom', 'nloom'});
addParameter(opt, 'type', 'cumulative', @(x) any(strcmp(x, {'cumulative', 'kde', 'whisk'})));
addParameter(opt, 'quantiles', 'group', @(x) any(strcmp(x, {'group', 'pooled'})));
addParameter(opt, 'n_quantiles', 3, @(x) isnumeric(x) && isscalar(x) && x >= 2 && x == round(x));
addParameter(opt, 'sloom_to_plot', 'all');
addParameter(opt, 'moving_to_plot', 'all');
addParameter(opt, 'check_quantiles', false);
parse(opt, varargin{:});

T               = opt.Results.bouts;
export          = opt.Results.export;
period          = opt.Results.period;
paths           = opt.Results.paths;
params          = opt.Results.params;
type            = opt.Results.type;
quantile_mode   = opt.Results.quantiles;
n_quantiles     = opt.Results.n_quantiles;
sloom_to_plot   = opt.Results.sloom_to_plot;
moving_to_plot  = opt.Results.moving_to_plot;
check_quantiles = opt.Results.check_quantiles;

if ischar(params) || isstring(params)
    params = cellstr(params);
end
if export && ~isstruct(paths)
    error('fd_distr_withparam_by_moving_flies:paths', 'export needs paths, a struct with the figure folder in .fig');
end

needed = {'durations_s', 'is_censored', 'ls', 'moving_flies', 'period'};
if ~istable(T) || isempty(T)
    error('fd_distr_withparam_by_moving_flies:noBouts', 'bouts has to be a table of freezes');
end
absent = needed(~ismember(needed, T.Properties.VariableNames));
if ~isempty(absent)
    error('fd_distr_withparam_by_moving_flies:columns', 'bouts has no column %s', strjoin(absent, ', '));
end

%% ===================== Infer period =====================
if isempty(period)
    u = unique(T.period);
    if numel(u) ~= 1
        error('Mixture of baseline and loom periods');
    end
    if u == 1
        period = 'loom';
    else
        period = 'bsl';
    end
end

%% ===================== Which freezes are plotted =====================
all_sloom = unique(T.ls)';
switch sloom_to_plot
    case 'all'
        sloom_vals = all_sloom;
    case 'slow'
        sloom_vals = min(all_sloom);
    case 'fast'
        sloom_vals = max(all_sloom);
    otherwise
        sloom_vals = sloom_to_plot;
end
T = T(ismember(T.ls, sloom_vals), :);

all_moving = unique(T.moving_flies)';
if ischar(moving_to_plot) && strcmp(moving_to_plot, 'all')
    levels = all_moving;
else
    levels = reshape(moving_to_plot, 1, []);
    levels = levels(ismember(levels, all_moving));
end
T = T(ismember(T.moving_flies, levels), :);
if isempty(T)
    error('fd_distr_withparam_by_moving_flies:noBouts', 'No freezes left for this loom speed and number of moving flies');
end

moving_cmap = cmapper([], 5).vars.moving_flies(2:end, :);

%% ===================== Loop over parameters =====================
for idx_param = params
    param = idx_param{1};
    if ~ismember(param, T.Properties.VariableNames)
        error('fd_distr_withparam_by_moving_flies:columns', 'bouts has no column %s', param);
    end
    values = T.(param);
    is_split = ~strcmp(param, 'moving_flies');

    [nq, thr_pooled, data_driven] = param_settings(param, T, n_quantiles);
    col  = cmapper([], nq);
    cmap = get_var_cmap(col, param);
    mid_q = ceil(nq / 2) + 1;

    if is_split
        tile_levels = levels;
    else
        tile_levels = NaN;
    end
    n_tiles = numel(tile_levels);

    %% Which level of the variable every freeze is in, from the edges of its panel
    tile_rows = cell(1, n_tiles);
    thr_tile = cell(1, n_tiles);
    q_all = nan(height(T), 1);
    for k = 1:n_tiles
        if is_split
            tile_rows{k} = T.moving_flies == tile_levels(k);
        else
            tile_rows{k} = true(height(T), 1);
        end
        thr_k = thr_pooled;
        if strcmp(quantile_mode, 'group') && data_driven
            thr_k = quantile(values(tile_rows{k}), linspace(0, 1, nq + 1));
        end
        if any(diff(thr_k) <= 0)
            error('fd_distr_withparam_by_moving_flies:edges', ...
                'The levels of %s collapse in the panel of %d moving flies (repeated values or too few freezes)', param, tile_levels(k));
        end
        thr_tile{k} = thr_k;
        q_all(tile_rows{k}) = discretize(values(tile_rows{k}), thr_k);
    end

    if check_quantiles
        fprintf('\n%s, %s quantiles\n', param, quantile_mode);
        for k = 1:n_tiles
            n_level = histcounts(q_all(tile_rows{k}), 0.5:1:nq + 0.5);
            if is_split, fprintf('  %d moving flies', tile_levels(k)); else, fprintf('  all freezes'); end
            fprintf(': edges [%s] | freezes per level [%s]\n', sprintf('%.4g ', thr_tile{k}), sprintf('%d ', n_level));
        end
    end

    %% Common bins and limits, so that the panels can be compared
    [bin_edges, x_lim] = control_bins(values, param);
    counts_max = 0;
    for k = 1:n_tiles
        if isempty(bin_edges)
            c = histcounts(q_all(tile_rows{k}), 0.5:1:nq + 0.5);
        else
            c = histcounts(values(tile_rows{k}), bin_edges);
        end
        counts_max = max([counts_max, c]);
    end

    %% ===================== Figure =====================
    if n_tiles == 1
        fig_width = 380;
    else
        fig_width = 150 + 190 * n_tiles;
    end
    fh = figure('Color', 'w', 'Position', [100 100 fig_width 460]);
    tl = tiledlayout(12, n_tiles, 'TileSpacing', 'tight', 'Padding', 'compact');

    %% ===================== Control distribution, one per panel =====================
    for k = 1:n_tiles
        ax = nexttile([3 1]); hold on
        plot_control_distribution(values(tile_rows{k}), q_all(tile_rows{k}), thr_tile{k}, cmap, param, bin_edges, x_lim);
        ylim([0 max(counts_max, 1)]);
        pad_ylim(ax, 0.075)
        apply_generic(ax, 'no_y', true, 'no_x', false, 'font_size', 20, 'tick_length', 0.03)
        xt = xticks(ax);
        if numel(xt) > 4
            xticks(ax, xt(1:2:end));
        end
        ax.XAxis.Exponent = 0;
    end

    %% ===================== Colours of the levels =====================
    name_label = control_label(param);
    if strcmp(quantile_mode, 'group') && data_driven && is_split
        name_label = [name_label ', ' quantile_name(nq) ' within each group'];
    end
    % The name has to fit in the width of the figure, which depends on the number of panels
    name_size = min(20, max(10, floor(0.9 * (fig_width - 30) / (numel(name_label) * 0.445))));
    nexttile([2 n_tiles]);
    draw_colour_strip(gca, cmap(2:end, :), quantile_labels(param, nq, T), name_label, cmap(mid_q, :), name_size);

    %% ===================== Duration distributions, one per panel =====================
    ax_bottom = gobjects(1, n_tiles);
    for k = 1:n_tiles
        ax = nexttile([7 1]); hold on
        plot_duration_distribution( ...
            T.durations_s(tile_rows{k}), q_all(tile_rows{k}), ...
            cmap, nq, type, period, logical(T.is_censored(tile_rows{k})));

        if k == 1
            ylabel(type);
            apply_generic(ax, 'tick_length', 0.025, 'font_size', 26, ...
                'yticks', [0 1], 'xticks', [0 10])
        else
            apply_generic(ax, 'tick_length', 0.025, 'font_size', 26, ...
                'xticks', [0 10], 'no_y', true)
        end
        if is_split
            add_group_marker(ax, tile_levels(k), moving_cmap(tile_levels(k) + 1, :), k == 1);
        end
        ax_bottom(k) = ax;
    end
    xlabel(tl, 'Freeze Duration (s)', 'FontSize', 26, 'FontName', 'Helvetica-Light');

    if numel(ax_bottom) > 1
        linkaxes(ax_bottom, 'y')
    end

    %% ===================== Export =====================
    if export
        exporter(fh, paths, sprintf('fd_%s_by_moving_flies_%s.pdf', param, period))
        exporter(fh, paths, sprintf('fd_%s_by_moving_flies_%s.png', param, period))
    end
end
end

%% =====================================================================
%% ============================ HELPERS ================================
%% =====================================================================

function [nq, thresholds, data_driven] = param_settings(param, T, n_quantiles)
% Number of levels and their edges. data_driven is true when the edges are quantiles
% of the data, false when they are fixed
data_driven = false;
switch param
    case {'nloom', 'ln'}
        % Edges between the loom numbers, so that 1-5, 6-10, 11-15 and 16-20 are the bins.
        % ln is nloom / 10
        nq = 4; thresholds = 0.5:5:20.5;
        if strcmp(param, 'ln'), thresholds = thresholds / 10; end
    case {'ls', 'jump_start', 'late_onset'}
        % binary: quantile() edges collapse to [0 0 0 0 1] and discretize
        % puts every bout in the top bin
        nq = 2; thresholds = -0.5:1:1.5;
    case 'moving_flies'
        nq = 5; thresholds = -0.5:1:5.5;
    otherwise
        nq = n_quantiles;
        thresholds = quantile(T.(param), linspace(0, 1, nq + 1));
        data_driven = true;
end
end

% ---------------------------------------------------------------------

function s = quantile_name(n)
names = {'halves', 'terciles', 'quartiles', 'quintiles'};
if n >= 2 && n <= 5
    s = names{n - 1};
else
    s = sprintf('%d quantiles', n);
end
end

% ---------------------------------------------------------------------

function s = control_label(param)
% What the control variable is called on the figure; the column name if it has no entry
switch param
    case 'avg_sm_freeze_norm', s = 'Avg. Social Motion';
    case 'avg_fs_1s_norm',     s = 'Focal Speed before Loom';
    case 'ls',                 s = 'Loom Speed';
    case 'moving_flies',       s = 'Moving Flies';
    case 'late_onset',         s = 'Onset after split';
    case {'nloom', 'ln'},      s = 'Loom Number';
    case 'sm_1s',              s = 'Social Motion (1st second)';
    case 'sm_2s',              s = 'Social Motion (1st 2 s)';
    case 'sm',                 s = 'Social Motion (freeze)';
    case 'avg_fs_1s',          s = 'Focal Speed (1 s before freeze)';
    case 'avg_fs_loom',        s = 'Focal Speed (before previous loom)';
    otherwise,                 s = strrep(param, '_', ' ');
end
end

% ---------------------------------------------------------------------

function [edges, x_lim] = control_bins(data, param)
% Histogram bins and x limits of a continuous variable, from all the panels together
edges = [];
x_lim = [];
% Discrete variables: no bins, one bar per level
if ismember(param, {'ln', 'ls', 'nloom', 'moving_flies', 'jump_start', 'late_onset'})
    return
end
data = data(~isnan(data));
iqr_val = iqr(data);
if iqr_val > 0
    bw = 2 * iqr_val / numel(data)^(1/2.1);
    nbins = max(20, ceil(range(data) / bw));
else
    nbins = 50;
end
[~, edges] = histcounts(data, nbins);

lo = prctile(data, 1);
hi = prctile(data, 99);
if lo < hi
    x_lim = [lo hi];
else
    x_lim = [min(data) max(data)];
end
end

% ---------------------------------------------------------------------

function plot_control_distribution(data, q, thresholds, cmap, param, bin_edges, x_lim)
% Distribution of the variable in one panel, each bin the colour of its level
if isempty(bin_edges)
    n_level = size(cmap, 1) - 1;
    bh = bar(1:n_level, histcounts(q, 0.5:1:n_level + 0.5), ...
        'FaceColor', 'flat', 'EdgeColor', 'none');
    bh.CData = cmap(2:end, :);
    xlim([.5 n_level + .5])
    xticks([])
else
    counts = histcounts(data, bin_edges);
    centers = bin_edges(1:end-1) + diff(bin_edges) / 2;

    b = bar(centers, counts, 1, 'FaceColor', 'flat', 'EdgeColor', 'none');
    b.CData = repmat(cmap(2, :), numel(centers), 1);
    for level = 2:numel(thresholds) - 1
        b.CData(centers >= thresholds(level), :) = repmat(cmap(1 + level, :), sum(centers >= thresholds(level)), 1);
    end
    xlim(x_lim)
end
end

% ---------------------------------------------------------------------

function plot_duration_distribution(durations, qmask, cmap, nq, type, period, censored)

switch type
    case 'cumulative'
        for q = 1:nq
            d = durations(qmask == q);
            c = censored(qmask == q);
            if isempty(d) || all(c), continue, end
            [f, x] = ecdf(d, 'Censoring', c);
            stairs(x, f, 'LineWidth', 3, 'Color', cmap(1 + q, :))
        end
        ylim([0 1])
        pad_ylim(gca, 0.025)

    case 'kde'
        % The Kaplan-Meier mass of every bin, so a freeze that was cut counts where it
        % was seen to last and not as if it ended there. A bin that reaches past the
        % longest observation, when that one was cut, is left blank, and a curve adds
        % up to the F(t) of the cumulative plot, not to 1
        bin_w = 6 / 60;
        n_bin = max(1, ceil((max(durations) - min(durations)) / bin_w));
        edges = min(durations) + (0:n_bin) * bin_w;
        for q = 1:nq
            d = durations(qmask == q);
            c = censored(qmask == q);
            if isempty(d) || all(c), continue, end
            [mass, covered] = km_bin_mass(d, c, edges);
            h = mass / bin_w;
            h(~covered) = NaN;
            stairs(edges, [h h(end)], 'LineWidth', 1.25, 'Color', cmap(1 + q, :))
        end
        pad_ylim(gca, 0.005)
        ylim([-0.05 1.05])

    case 'whisk'
        % Does not take the censoring into account
        for q = 1:nq
            d = durations(qmask == q);
            if isempty(d), continue, end
            boxchart(q * ones(size(d)), d, ...
                'Orientation', 'horizontal', ...
                'BoxWidth', 0.5, ...
                'BoxFaceColor', cmap(1 + q, :), ...
                'MarkerStyle', 'none', ...
                'LineWidth', 1.2);
        end
        yticks(1:nq)
        yticklabels(arrayfun(@(x) ['Q' num2str(x)], 1:nq, 'uni', 0))
        ylim([0.5 nq + 0.5])
end

switch period
    case 'loom', xlim([0 10])
    case 'bsl',  xlim([-0.05 1.05])
end
end

% ---------------------------------------------------------------------

function [mass, covered] = km_bin_mass(d, c, edges)
% Kaplan-Meier mass of every bin: the steps of the curve that fall in it. covered is
% false for the bins that reach past the longest observation when that one was cut,
% their mass is not known
covered = true(1, numel(edges) - 1);
[f, x] = ecdf(d, 'Censoring', c);
step = diff(f);
b = discretize(x(2:end), edges);
ok = ~isnan(b);
mass = accumarray(b(ok), step(ok), [numel(edges) - 1, 1])';
if f(end) < 1 - 1e-9
    covered = edges(2:end) <= max(d) + 1e-9;
end
end

% ---------------------------------------------------------------------

function cmap = get_var_cmap(col, varname)
if isfield(col.vars, varname)
    cmap = col.vars.(varname);
else
    cmap = col.vars.default;
end
end

% ---------------------------------------------------------------------

function labels = quantile_labels(param, nq, T)
dash = char(8211);
switch param
    case {'nloom', 'ln'}
        % The first bin starts at the first loom in the data, 2 when loom 1 is left out
        labels = {sprintf('%d%c5', max(1, min(T.nloom)), dash), ['6' dash '10'], ['11' dash '15'], ['16' dash '20']};
    case {'sloom', 'ls'}
        labels = {'Slow', 'Fast'};
    case 'late_onset'
        labels = {'Early', 'Late'};
    case 'moving_flies'
        labels = arrayfun(@num2str, 0:nq - 1, 'uni', 0);
    otherwise
        labels = arrayfun(@(x) ['Q' num2str(x)], 1:nq, 'uni', 0);
end
end

% ---------------------------------------------------------------------

function draw_colour_strip(ax, colours, labels, name_label, name_colour, name_size)
% The colours of the levels in a bar with their labels below and the name of the
% variable above, across all the panels
nq = size(colours, 1);
hold(ax, 'on')
% image takes the centres of the first and last pixels: two identical rows from 0 to 1
% and one column per level from 0.5 to nq + 0.5
image(ax, [1 nq], [0.25 0.75], repmat(reshape(colours, 1, nq, 3), 2, 1));
ax.YDir = 'normal';
rectangle(ax, 'Position', [0.5 0 nq 1], 'EdgeColor', 'k', 'LineWidth', 2);
for k = 1:nq
    text(ax, k, -0.5, labels{k}, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', 'FontSize', 22);
end
text(ax, (nq + 1) / 2, 1.6, name_label, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
    'FontSize', name_size, 'Color', name_colour, 'Interpreter', 'none');
xlim(ax, [0.5 nq + 0.5]);
ylim(ax, [-3 4]);
axis(ax, 'off')
end

% ---------------------------------------------------------------------

function add_group_marker(ax, n_moving, face_colour, first)
% The number of moving flies of a panel, in a circle of its colour, in data coordinates
xl = xlim(ax);
yl = ylim(ax);
x = xl(1) + 0.125 * diff(xl);
y = yl(2) - 0.10 * diff(yl);
scatter(ax, x, y, 450, 'filled', 'MarkerFaceColor', face_colour, 'MarkerEdgeColor', 'none', 'HandleVisibility', 'off');
text(ax, x, y, sprintf('%d', n_moving), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'FontSize', 18);
if first
    text(ax, x + 0.12 * diff(xl), y, 'moving flies', 'FontSize', 14, 'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle');
end
end

% ---------------------------------------------------------------------

function pad_ylim(ax, fraction)
yl = ax.YLim;
p = fraction * diff(yl);
ax.YLim = [yl(1) - p yl(2) + p];
end
