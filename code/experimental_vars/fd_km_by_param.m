function [fh, stats] = fd_km_by_param(varargin)
% FD_KM_BY_PARAM  Kaplan-Meier freeze-exit curves stratified by a baseline
% covariate (default: social motion), with Greenwood confidence bands, a
% number-at-risk table and a Cox hazard ratio per panel.
%
%   [fh, stats] = fd_km_by_param('bouts', T, 'param', 'sm_early')
%
% T comes from build_freeze_survival: it needs t_obs (frames), event, ls and
% the covariate column, and its UserData carries the landmark. Every bout in T
% is at risk from the landmark on, so the curves are
%
%       F(t) = P(freeze has ended by t | freeze lasted past the landmark,
%                                        no contact before it)
%
% with contact right-censored. They are the contact-free freeze duration and
% not the observed one, see the caveat about dependent censoring in the notes
% of build_freeze_survival. The covariate must be measured before the
% landmark, otherwise it is partly the outcome.
%
% Options
%   bouts     table
%   param     'sm_early'  covariate column
%   nq        4           quantile strata, edges taken on the pooled data so
%                         the strata mean the same thing in every panel
%   split_by  'ls'        one panel per value; '' for a single panel
%   n_min     30          the time axis stops where the thinnest stratum has
%                         fewer bouts than this at risk
%   xmax_s    10          upper limit of the time axis (s)
%   ci        true        Greenwood bands, log-log, 95%
%   export    false       with 'paths' (path_generator output) and 'name'
%
% stats (struct array, one element per panel): n, n_events, n_contact,
% hr_per_sd, hr_lo, hr_hi, z, xmax_s, at_risk (strata x ticks), ticks_s.
%
% The Cox hazard ratio is per +1 SD of the continuous covariate, within the
% panel. It assumes proportional hazards over the strata, which the curves
% let you eyeball.

opt = inputParser;
addParameter(opt, 'bouts', []);
addParameter(opt, 'param', 'sm_early');
addParameter(opt, 'nq', 4);
addParameter(opt, 'split_by', 'ls');
addParameter(opt, 'n_min', 30);
addParameter(opt, 'xmax_s', 10);
addParameter(opt, 'ci', true);
addParameter(opt, 'export', false);
addParameter(opt, 'paths', []);
addParameter(opt, 'name', '');
addParameter(opt, 'font_size', 16);
parse(opt, varargin{:});

T     = opt.Results.bouts;
param = opt.Results.param;
nq    = opt.Results.nq;
n_min = opt.Results.n_min;
fs    = opt.Results.font_size;
z_ci  = 1.96;

ud = T.Properties.UserData;
if isstruct(ud) && isfield(ud, 'landmark'), origin_f = ud.landmark; else, origin_f = 30; end
origin_s = origin_f / 60;

T = T(~isnan(T.(param)), :);
x = T.(param);
edges = quantile(x, linspace(0, 1, nq + 1));
stratum = discretize(x, edges);

% Colours: social-motion variables share one palette
col = cmapper([], nq);
if startsWith(param, 'sm') && isfield(col.vars, 'sm')
    cmap = col.vars.sm;
elseif isfield(col.vars, param)
    cmap = col.vars.(param);
else
    cmap = col.vars.default;
end
cmap = cmap(2:end, :);

if isempty(opt.Results.split_by)
    panel_vals = 1;  panel_id = ones(height(T), 1);
else
    panel_id = T.(opt.Results.split_by);
    panel_vals = unique(panel_id)';
end
nP = numel(panel_vals);

%% Layout: distribution of the covariate, curves, at-risk table
fh = figure('Color', 'w', 'Position', [100 100 330 * nP + 120, 700]);
tiledlayout(13, nP, 'TileSpacing', 'compact', 'Padding', 'compact');

ax_d = nexttile(1, [3 nP]); hold on
bw = 2 * iqr(x) / numel(x)^(1/3);
[c, e] = histcounts(x, max(20, ceil(range(x) / bw)));
ctr = e(1:end-1) + diff(e) / 2;
bh = bar(ctr, c, 1, 'FaceColor', 'flat', 'EdgeColor', 'none');
q_of_bin = discretize(ctr, edges);
q_of_bin(isnan(q_of_bin)) = nq;
bh.CData = cmap(q_of_bin, :);
xlim([prctile(x, 0.5) prctile(x, 99.5)])
apply_generic(ax_d, 'no_y', true, 'font_size', fs)
if strcmp(param, 'sm_early')
    xlabel(sprintf('Social motion, first %g s of the freeze', origin_s))
else
    xlabel(strrep(param, '_', ' '), 'Interpreter', 'none')
end

stats = repmat(struct('n', [], 'n_events', [], 'n_contact', [], 'hr_per_sd', [], 'hr_lo', [], ...
    'hr_hi', [], 'z', [], 'xmax_s', [], 'at_risk', [], 'ticks_s', []), 1, nP);
ax_c = gobjects(1, nP);

for ip = 1:nP
    in_panel = panel_id == panel_vals(ip);
    t  = T.t_obs(in_panel);
    ev = T.event(in_panel);
    st = stratum(in_panel);

    % Time axis ends where the thinnest stratum runs out of bouts at risk
    grid_f = origin_f : opt.Results.xmax_s * 60;
    at_risk_grid = zeros(numel(grid_f), nq);
    for q = 1:nq
        at_risk_grid(:, q) = arrayfun(@(g) sum(t(st == q) >= g), grid_f);
    end
    last_ok = find(min(at_risk_grid, [], 2) >= n_min, 1, 'last');
    if isempty(last_ok), last_ok = 1; end
    xmax_f = grid_f(last_ok);

    ax_c(ip) = nexttile(3 * nP + ip, [7 1]); hold on
    for q = 1:nq
        m = st == q;
        if ~any(m), continue, end
        km = km_fit(t(m), ev(m), z_ci);
        keep = km.t <= xmax_f;
        xs = [origin_f; km.t(keep)] / 60;
        F  = [0; 1 - km.S(keep)];
        if opt.Results.ci
            Flo = [0; 1 - km.hi(keep)];
            Fhi = [0; 1 - km.lo(keep)];
            ok = ~isnan(Flo) & ~isnan(Fhi);
            [xb, ylo] = stairs(xs(ok), Flo(ok));
            [~,  yhi] = stairs(xs(ok), Fhi(ok));
            fill([xb; flipud(xb)], [ylo; flipud(yhi)], cmap(q, :), ...
                'FaceAlpha', 0.18, 'EdgeColor', 'none');
        end
        stairs(xs, F, 'LineWidth', 3, 'Color', cmap(q, :));
    end
    ylim([-0.025 1.025])
    xlim([origin_s xmax_f / 60])
    apply_generic(ax_c(ip), 'font_size', fs, 'tick_length', 0.02, 'line_width', 2, 'yticks', [0 0.5 1])
    if ip == 1, ylabel('P(freeze ended)'); else, set(ax_c(ip), 'YTickLabel', []); end
    set(ax_c(ip), 'XTickLabel', [])

    if ~isempty(opt.Results.split_by)
        lbl = sprintf('%s = %g', opt.Results.split_by, panel_vals(ip));
        if strcmp(opt.Results.split_by, 'ls')
            lbl = ternary(panel_vals(ip) == min(panel_vals), 'Slow loom', 'Fast loom');
        end
        title(lbl, 'FontWeight', 'normal', 'FontSize', fs + 2)
    end

    % Cox hazard ratio per SD of the continuous covariate
    xz = zscore(T.(param)(in_panel));
    [b, ~, ~, s] = coxphfit(xz, t, 'Censoring', double(~ev));
    hr = exp(b); lo = exp(b - z_ci * s.se); hi = exp(b + z_ci * s.se);
    text(ax_c(ip), xmax_f / 60, 0.03, sprintf('HR %.2f [%.2f, %.2f] per SD', hr, lo, hi), ...
        'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom', 'FontSize', fs - 2);

    % At-risk table
    ticks_s = unique([origin_s, 2:2:floor(xmax_f / 60)]);
    ticks_s = ticks_s(ticks_s >= origin_s);
    ar = zeros(nq, numel(ticks_s));
    for q = 1:nq
        for k = 1:numel(ticks_s)
            ar(q, k) = sum(t(st == q) >= ticks_s(k) * 60);
        end
    end
    set(ax_c(ip), 'XTick', ticks_s)
    set(ax_c(ip), 'XTickLabel', arrayfun(@(v) sprintf('%g', v), ticks_s, 'UniformOutput', false))
    xlabel(ax_c(ip), 'Time since freeze onset (s)')

    ax_r = nexttile(10 * nP + ip, [3 1]); hold on
    xlim(ax_r, [origin_s xmax_f / 60]); ylim(ax_r, [0.5 nq + 1.3]); axis(ax_r, 'off')
    ax_r.Clipping = 'off';
    for q = 1:nq
        for k = 1:numel(ticks_s)
            text(ax_r, ticks_s(k), nq + 1 - q, sprintf('%d', ar(q, k)), 'Color', cmap(q, :), ...
                'HorizontalAlignment', 'center', 'FontSize', fs - 4, 'FontWeight', 'bold');
        end
        if ip == 1
            text(ax_r, origin_s - 0.07 * (xmax_f / 60 - origin_s), nq + 1 - q, sprintf('Q%d', q), ...
                'Color', cmap(q, :), 'HorizontalAlignment', 'right', 'FontSize', fs - 4, 'FontWeight', 'bold');
        end
    end
    if ip == 1
        text(ax_r, origin_s - 0.07 * (xmax_f / 60 - origin_s), nq + 1.2, 'No. at risk', ...
            'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom', 'FontSize', fs - 4);
    end

    stats(ip) = struct('n', sum(in_panel), 'n_events', sum(ev), ...
        'n_contact', sum(T.censored_contacts(in_panel)), 'hr_per_sd', hr, 'hr_lo', lo, ...
        'hr_hi', hi, 'z', b / s.se, 'xmax_s', xmax_f / 60, 'at_risk', ar, 'ticks_s', ticks_s);
end
linkaxes(ax_c, 'y')

if opt.Results.export
    exporter(fh, opt.Results.paths, opt.Results.name);
end
end

% ---------------------------------------------------------------------

function km = km_fit(t, ev, z)
% Product-limit estimate with Greenwood variance and log-log bands.
% Events precede censorings at tied times (the risk set at t is #{t_j >= t}),
% the usual convention and the one ecdf(..., 'Censoring') uses
t = t(:); ev = logical(ev(:));
ut = unique(t(ev));
nr = arrayfun(@(u) sum(t >= u), ut);
d  = arrayfun(@(u) sum(t == u & ev), ut);
S  = cumprod(1 - d ./ nr);
G  = cumsum(d ./ (nr .* (nr - d)));
se = sqrt(G) ./ abs(log(S));
se(~isfinite(se)) = NaN;
km.t  = ut;
km.S  = S;
km.lo = S .^ exp(z * se);
km.hi = S .^ exp(-z * se);
end

function out = ternary(cond, a, b)
if cond, out = a; else, out = b; end
end
