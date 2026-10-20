function [fh, stats] = fd_km_grid(varargin)
% FD_KM_GRID  Grid of freeze-duration distributions with the censoring taken into
% account, with any variables on the rows and on the columns.
%
%   [fh, stats] = fd_km_grid('bouts', freezes)
%   [fh, stats] = fd_km_grid('bouts', freezes, 'rows', {'sm_1s'}, 'cols', {'avg_fs_1s', 'ls'})
%   [fh, stats] = fd_km_grid('bouts', freezes, 'rows', {{'avg_fs_1s', 4}}, 'cols', {'ls'})
%
% Each panel is a histogram of the freeze duration made of the Kaplan-Meier mass
% of every bin, as the 'kde' type of fd_distr_withparam_new, so a freeze that was
% cut by a contact or by the end of the loom cycle counts where it was seen to
% last and not as if it had ended there. quantilizer_v2 and ksdensity do not
% know about the censoring: they count the freezes that were cut as ended ones,
% which piles up mass at the time of the cut and shortens the tail.
%
% bouts is the table of loom_evoked_freezes, with durations_s, is_censored and the
% variables asked for. Bins beyond the longest observation are left blank when
% that observation was cut, so the bars of a panel add up to P(freeze ended by the
% end of the axis), not to 1.
%
% Options
%   rows, cols  {'sm_1s'}, {'avg_fs_1s', 'ls'}   cell arrays of variables. A grid
%                 line is a combination of their levels, the last variable
%                 changing fastest, so the default columns are slow / fast loom
%                 at the low focal speed, then at the high one. A variable is a
%                 column name, or {name, n_levels} for n quantile levels of a
%                 continuous one (default 3 for the first variable of the rows,
%                 2 for any other). Fixed levels: ls (slow, fast), nloom
%                 (2-5, 6-10, 11-15, 16-20), moving_flies, late_onset, jump_start
%   edges       'nested'   where the quantiles of a continuous variable are taken:
%                 'nested' within every combination of the variables before it
%                          (the fixed ones first, then the continuous ones in
%                          the order cols, rows), so the cells of a line have
%                          equal numbers and every variable is compared at the
%                          same level of the others
%                 'cat'    within the fixed variables only, so continuous
%                          variables are not conditioned on each other
%                 'pooled' on all the freezes (quantilizer_v2 'global')
%   type        'hist'     'hist' bars, 'cumulative' Kaplan-Meier curves
%   bin_f       10         bin width in frames
%   xmax_s      10         upper limit of the time axis
%   min_n       20         a panel with fewer freezes is left empty
%   export      false      with 'paths' and 'name'
%
% Colours are one per row. stats(row, col) has n, n_cens, frac_cens, s_end (S at
% xmax_s, NaN if the follow-up does not reach it), levels (the level of every
% variable, in the order cols then rows) and range (the min and max of every
% continuous variable in the cell).

opt = inputParser;
addParameter(opt, 'bouts', []);
addParameter(opt, 'rows', {'sm_1s'});
addParameter(opt, 'cols', {'avg_fs_1s', 'ls'});
addParameter(opt, 'edges', 'nested', @(x) any(strcmp(x, {'nested', 'cat', 'pooled'})));
addParameter(opt, 'type', 'hist', @(x) any(strcmp(x, {'hist', 'cumulative'})));
addParameter(opt, 'bin_f', 10);
addParameter(opt, 'xmax_s', 10);
addParameter(opt, 'min_n', 20);
addParameter(opt, 'font_size', 18);
addParameter(opt, 'export', false);
addParameter(opt, 'paths', []);
addParameter(opt, 'name', 'km_grid.pdf');
parse(opt, varargin{:});

T     = opt.Results.bouts;
xmax  = opt.Results.xmax_s;
bin_w = opt.Results.bin_f / 60;
fs    = opt.Results.font_size;
col   = cmapper();

row_spec = parse_vars(opt.Results.rows, 3);
col_spec = parse_vars(opt.Results.cols, 2);
specs = [col_spec, row_spec];
nv = numel(specs);
nR = prod([row_spec.n]);
nC = prod([col_spec.n]);

T = T(~any(isnan(cell2mat(cellfun(@(v) T.(v), {specs.name}, 'UniformOutput', false))), 2), :);

%% Level of every freeze on every variable
level = nan(height(T), nv);
fixed = [specs.cat];
for v = find(fixed)
    level(:, v) = fixed_level(T, specs(v));
end
order = [find(fixed), find(~fixed)];
assigned = [];
for v = order
    if fixed(v), assigned(end + 1) = v; continue, end %#ok<AGROW>
    switch opt.Results.edges
        case 'pooled', ctx = [];
        case 'cat',    ctx = assigned(fixed(assigned));
        case 'nested', ctx = assigned;
    end
    if isempty(ctx), key = ones(height(T), 1); else, [~, ~, key] = unique(level(:, ctx), 'rows'); end
    x = T.(specs(v).name);
    for k = unique(key)'
        in_k = key == k;
        e = quantile(x(in_k), linspace(0, 1, specs(v).n + 1));
        level(in_k, v) = discretize(x(in_k), e);
    end
    assigned(end + 1) = v; %#ok<AGROW>
end

% Grid position, the last variable of a list changing fastest
ncv = numel(col_spec);
row_idx = combo_index(level(:, ncv + 1:end), [row_spec.n]);
col_idx = combo_index(level(:, 1:ncv), [col_spec.n]);
col_titles = combo_labels(col_spec);
row_titles = combo_labels(row_spec);

%% Figure
fh = figure('Color', 'w', 'Position', [100 100 280 * nC + 100, 230 * nR + 100]);
tl = tiledlayout(nR, nC, 'TileSpacing', 'compact', 'Padding', 'compact');
ax = gobjects(nR, nC);
stats = repmat(struct('n', 0, 'n_cens', 0, 'frac_cens', NaN, 's_end', NaN, 'levels', [], 'range', []), nR, nC);
edges_t = 0:bin_w:xmax;
ymax = 0;
row_col = row_colors(col, nR);

for r = 1:nR
    for j = 1:nC
        ax(r, j) = nexttile(tl, (r - 1) * nC + j); hold on
        m = row_idx == r & col_idx == j;
        d = T.durations_s(m);
        c = logical(T.is_censored(m));
        cc = row_col(r, :);
        stats(r, j).n = numel(d);
        stats(r, j).n_cens = sum(c);
        stats(r, j).frac_cens = mean(c);
        stats(r, j).levels = level(find(m, 1), :);
        stats(r, j).range = cell2mat(arrayfun(@(v) [min(T.(specs(v).name)(m)) max(T.(specs(v).name)(m))], ...
            find(~fixed), 'UniformOutput', false));

        if numel(d) >= opt.Results.min_n && ~all(c)
            [f, x] = ecdf(d, 'Censoring', c);
            if max(d) >= xmax || f(end) >= 1 - 1e-9
                stats(r, j).s_end = 1 - interp1(x(2:end), f(2:end), xmax, 'previous', f(end));
            end
            switch opt.Results.type
                case 'hist'
                    [mass, covered] = km_bin_mass(d, c, edges_t);
                    h = mass / bin_w;
                    h(~covered) = NaN;
                    centers = edges_t(1:end-1) + bin_w / 2;
                    bar(centers(covered), h(covered), 1, 'FaceColor', cc, 'FaceAlpha', 0.35, 'EdgeColor', 'none', 'ShowBaseLine', 'off');
                    stairs(edges_t, [h h(end)], 'LineWidth', 2, 'Color', cc);
                    ymax = max([ymax, h]);
                case 'cumulative'
                    stairs(x, f, 'LineWidth', 3, 'Color', cc);
            end
        end
        if r == 1, title(col_titles{j}, 'FontWeight', 'normal', 'FontSize', fs - 2); end
        xlim([0 xmax])
        apply_generic(ax(r, j), 'font_size', fs, 'tick_length', 0.025, 'line_width', 2.5, 'xticks', [0 xmax / 2 xmax]);
        if j > 1, set(ax(r, j), 'YTickLabel', []); ax(r,j).YAxis.Visible =  'off'; end

    end
end

if strcmp(opt.Results.type, 'hist')
    ylim(ax(:), [-0.025 1.05 * ymax]);  ylab = 'KM mass / s';
else
    ylim(ax(:), [0 1]);            ylab = 'P(freeze ended)';
end
linkaxes(ax(:), 'xy');
xlabel(tl, 'Freeze duration (s)', 'FontSize', fs + 2);
ylabel(tl, ylab, 'FontSize', fs + 2);

% Counts, in data coordinates once the limits are set
yl = ylim(ax(1));
for r = 1:nR
    for j = 1:nC
        text(ax(r, j), 0.98 * xmax, yl(1) + 0.95 * diff(yl), ...
            sprintf('n=%d (%.0f%% cens.)', stats(r, j).n, 100 * stats(r, j).frac_cens), ...
            'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', 'FontSize', fs - 6);
    end
    ylabel(ax(r, 1), row_titles{r}, 'FontSize', fs - 2);
end

if opt.Results.export
    exporter(fh, opt.Results.paths, opt.Results.name);
end
end

% ---------------------------------------------------------------------

function s = parse_vars(list, n_first)
% {'a', {'b', 4}} to a struct array of name, n, cat. The first variable of the list
% has n_first levels by default, the others 2
fixed_vars = {'ls', 'nloom', 'moving_flies', 'late_onset', 'jump_start'};
if ~iscell(list), list = {list}; end
for k = 1:numel(list)
    v = list{k};
    if iscell(v), name = v{1}; n = v{2}; else, name = v; n = ternary(k == 1, n_first, 2); end
    s(k).name = name; %#ok<AGROW>
    s(k).cat = ismember(name, fixed_vars);
    s(k).n = n;
    if s(k).cat
        s(k).n = fixed_n(name);
    end
end
end

function n = fixed_n(name)
switch name
    case 'nloom', n = 4;
    case 'moving_flies', n = 5;
    otherwise, n = 2;
end
end

function lv = fixed_level(T, spec)
x = T.(spec.name);
switch spec.name
    case 'nloom',        lv = discretize(x, 0.5:5:20.5);
    case 'moving_flies', lv = discretize(x, -0.5:1:5.5);
    otherwise,           lv = discretize(x, [-0.5 0.5 1.5]);   % ls 0 slow, 1 fast; 0/1 flags
end
end

function idx = combo_index(lv, ns)
% Linear index of the combination of levels, the last column changing fastest
idx = ones(size(lv, 1), 1);
for k = 1:numel(ns)
    idx = (idx - 1) * ns(k) + lv(:, k);
end
end

function labels = combo_labels(spec)
% One label per combination of levels of the variables, last one fastest
n = [spec.n];
labels = cell(1, prod(n));
for i = 1:prod(n)
    lv = zeros(1, numel(n));
    rem_i = i - 1;
    for k = numel(n):-1:1
        lv(k) = mod(rem_i, n(k)) + 1;
        rem_i = floor(rem_i / n(k));
    end
    parts = arrayfun(@(k) level_label(spec(k), lv(k)), 1:numel(n), 'UniformOutput', false);
    labels{i} = strjoin(parts, newline);
end
end

function s = level_label(spec, l)
switch spec.name
    case 'ls',           s = ternary(l == 1, 'Slow loom', 'Fast loom');
    case 'late_onset',   s = ternary(l == 1, 'Early', 'Late');
    case 'jump_start',   s = ternary(l == 1, 'No jump', 'Jump');
    case 'nloom',        names = {'Looms 2-5', 'Looms 6-10', 'Looms 11-15', 'Looms 16-20'}; s = names{l};
    case 'moving_flies', s = sprintf('%d moving', l - 1);
    otherwise
        short = struct('sm_1s', 'SM', 'sm_2s', 'SM 2 s', 'sm', 'SM', 'avg_fs_1s', 'Speed', 'avg_fs_loom', 'Speed');
        nm = ternary(isfield(short, spec.name), @() short.(spec.name), @() strrep(spec.name, '_', ' '));
        if spec.n == 2, names = {'low', 'high'}; s = [nm ' ' names{l}];
        elseif spec.n == 3, names = {'low', 'mid', 'high'}; s = [nm ' ' names{l}];
        else, s = sprintf('%s Q%d', nm, l);
        end
end
end

function cm = row_colors(col, nR)
base = col.empirical.col(1:2:end, :);   % one colour per row, as the first focal-speed half of the old grid
if nR <= size(base, 1), cm = base(1:nR, :); else, cm = lines(nR); end
end

function out = ternary(cond, a, b)
if cond, out = a; else, out = b; end
if isa(out, 'function_handle'), out = out(); end
end

function [mass, covered] = km_bin_mass(d, c, edges)
% Kaplan-Meier mass of every bin: the steps of the curve that fall in it. covered
% is false for the bins that reach past the longest observation when that one was
% cut, their mass is not known (same as in fd_distr_withparam_new)
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
