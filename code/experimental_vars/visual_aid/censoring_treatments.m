% Three ways of treating a right-censored freeze, on simulated data.
%
%   as recorded   the freeze counts as ending at the time it was cut, which is what
%                 the durations column holds
%   omit          the censored freezes are dropped
%   Kaplan-Meier  ecdf(..., 'Censoring', ...), what fd_distr_withparam_new and
%                 fd_distr_withparam_by_moving_flies plot
%
% The true durations are drawn here, so every estimate is compared with the
% empirical CDF of the full sample, which is what it tries to recover. The same n
% freezes are cut in three ways:
%
%   A  one freeze that lasted past 5 s is cut at 5 s
%   B  every freeze that lasts past 5 s is cut at 5 s, like the 10.5 s cap of
%      impose_contact_threshold
%   C  every freeze has its own cut time, drawn at random, like a collision
%
% Kaplan-Meier. With t_1 < t_2 < ... the times at which a freeze was seen to end,
%
%       n_i   freezes still under observation just before t_i (t_obs >= t_i, the
%             ones cut later included)
%       d_i   freezes that ended at t_i
%
%       S(t) = prod over t_i <= t of (1 - d_i / n_i)         F(t) = 1 - S(t)
%
% A cut freeze is in n_i until its cut time and in no step after it: it says "still
% frozen up to here" and nothing more. A freeze cut at the time another one ends is
% still in n_i for that step. The product is worked out step by step below
% (km_by_hand) and checked against ecdf.
%
% Figure: the CDFs on top, and below each the estimate minus the CDF of the full
% sample (the vertical scales of the bottom row are not shared). Ticks under the top
% panels mark the cut freezes. Kaplan-Meier stops at the longest observation when
% that one was cut, there is no information beyond it.
%
% The numbers printed at the end are the mean error of each estimate at three times,
% over n_rep samples, against the exact CDF of the exponential the durations are
% drawn from.
%
% Histograms. The same three treatments as the mass of every bin of width bin_w, in a
% second figure with the same layout. Kaplan-Meier has a mass at every step of its
% curve,
%
%       dF_i = S(t_i-) * d_i / n_i        S(t_i-) the survival just before t_i
%
% and the mass of a bin is the sum of the steps that fall in it. A cut freeze has no
% mass of its own: its 1/n goes to the freezes seen to end after its cut time, each
% of which weighs S / n_i instead of 1/n. A bin that reaches past the longest
% observation, when that one was cut, has an unknown mass and is left blank, so the
% histogram adds up to less than 1 there. The densities are mass / bin_w, in 1/s.
% The last section prints the same checks as above for the bin masses.

rng(4)
n      = 1000;      % freezes
mu     = 2;        % s, mean of the true duration (exponential)
t_cap  = 5;        % s, where the freezes are cut in A and B
mu_c   = 6;        % s, mean of the random cut time in C
n_rep  = 500;      % samples of the error table
t_grid = linspace(0, 10, 1001);
bin_w  = 0.25;      % s, bin width of the histograms
edges  = 0:bin_w:10;

scenario_name = {'A. one freeze cut at 5 s', 'B. all freezes past 5 s cut at 5 s', 'C. random cut times'};
method_name   = {'as recorded', 'omit', 'Kaplan-Meier'};
method_col    = [0.90 0.50 0.10; 0.20 0.50 0.80; 0 0 0];
line_width    = [1.5 1.5 1.5];
truth_col     = [0.75 0.75 0.75];

%% Figure
fh = figure('Color', 'w', 'Position', [100 100 900 450]);
tl = tiledlayout(2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

for s = 1:3
    [t_obs, cens, T] = simulate(s, n, mu, t_cap, mu_c);
    [f_true, x_true] = ecdf(T);
    F_true = eval_cdf(x_true, f_true, t_grid);

    ax_top = nexttile(s); hold on
    stairs(x_true, f_true, 'Color', truth_col, 'LineWidth', 7, 'DisplayName', 'full sample, no cut');
    ax_dif = nexttile(3 + s); hold on
    yline(ax_dif, 0, 'Color', truth_col, 'LineWidth', 2, 'HandleVisibility', 'off');

    % omit first and thicker, since it lies under as recorded where the two agree
    for m = [2 1 3]
        [x, f, last] = estimate(m, t_obs, cens);
        if isfinite(last)
            % The curve is flat from the last step to the longest observation
            stairs(ax_top, [x; last], [f; f(end)], 'Color', method_col(m, :), 'LineWidth', line_width(m), 'DisplayName', method_name{m});
        else
            stairs(ax_top, x, f, 'Color', method_col(m, :), 'LineWidth', line_width(m), 'DisplayName', method_name{m});
        end

        F = eval_cdf(x, f, t_grid);
        F(t_grid > last) = NaN;
        plot(ax_dif, t_grid, F - F_true, 'Color', method_col(m, :), 'LineWidth', line_width(m));
    end
    yl = ylim(ax_dif);
    ylim(ax_dif, [-.1 .1]);

    plot(ax_top, t_obs(cens), -0.04 * ones(nnz(cens), 1), '|', 'Color', [0.8 0.1 0.1], 'MarkerSize', 10, 'HandleVisibility', 'off');
    title(ax_top, {scenario_name{s}, sprintf('%d of %d cut', nnz(cens), n)}, 'FontWeight', 'normal');
    xlim(ax_top, [0 10]); ylim(ax_top, [-0.08 1.02]);
    xlim(ax_dif, [0 10]);
    apply_generic(ax_top, 'font_size', 16, 'tick_length', 0.025, 'yticks', [0 1], 'xticks', [0 5 10], 'line_width', 2);
    apply_generic(ax_dif, 'font_size', 16, 'tick_length', 0.025, 'xticks', [0 5 10], 'line_width', 2);
    if s == 1
        ylabel(ax_top, 'F(t)');
        ylabel(ax_dif, 'estimate - full sample');
        legend(ax_top, 'Location', 'southeast', 'Box', 'off', 'FontSize', 13);
    end
end
xlabel(tl, 'Freeze duration (s)', 'FontSize', 20);

%% The product worked out, around the cut freeze of A
rng(1)
[t_obs, cens] = simulate(1, n, mu, t_cap, mu_c);
tab = km_by_hand(t_obs, cens);
k = find(tab.c > 0);
fprintf('\nA, Kaplan-Meier by hand around the freeze cut at %g s (c = number cut at that time)\n', t_cap);
disp(tab(max(1, k - 3):min(height(tab), k + 3), :))
fprintf('step before the cut: 1/n = %.4f.   steps after: S(t_cut) / n_i = %.4f / %d = %.4f\n', ...
    1 / n, tab.S(k), tab.n_risk(k + 1), tab.S(k) / tab.n_risk(k + 1));

%% By hand against ecdf
rng(1)
fprintf('\nlargest difference between the product above and ecdf(..., ''Censoring'', ...)\n');
for s = 1:3
    [t_obs, cens] = simulate(s, n, mu, t_cap, mu_c);
    tab = km_by_hand(t_obs, cens);
    [x, f] = estimate(3, t_obs, cens);
    fprintf('  %-40s %.2e\n', scenario_name{s}, max(abs(tab.F - eval_cdf(x, f, tab.t))));
end

%% Error of each estimate over many samples
t_eval = [2 4 6];
F_exact = 1 - exp(-t_eval / mu);
est = nan(n_rep, 3, numel(t_eval));
fprintf('\nmean of (estimate - exact CDF) over %d samples of %d freezes\n', n_rep, n);
fprintf('exact CDF                                  %s\n', sprintf('%9.3f', F_exact));
fprintf('%-40s %s\n', '', sprintf('   t = %d s', t_eval));
for s = 1:3
    for r = 1:n_rep
        [t_obs, cens] = simulate(s, n, mu, t_cap, mu_c);
        for m = 1:3
            [x, f, last] = estimate(m, t_obs, cens);
            F = eval_cdf(x, f, t_eval);
            F(t_eval > last) = NaN;
            est(r, m, :) = F;
        end
    end
    bias = squeeze(mean(est, 1, 'omitnan')) - F_exact;
    fprintf('%s\n', scenario_name{s});
    for m = 1:3
        fprintf('  %-38s', method_name{m});
        for j = 1:numel(t_eval)
            if isnan(bias(m, j)), fprintf('%9s', 'n/a'); else, fprintf('%+9.3f', bias(m, j)); end
        end
        fprintf('\n');
    end
end
fprintf('n/a: the longest observation was cut before that time, Kaplan-Meier is not defined there\n');

%% Histogram figure
fh = figure('Color', 'w', 'Position', [100 100 900 450]);
tl = tiledlayout(2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
ax_top = gobjects(1, 3);
cut_t  = cell(1, 3);

for s = 1:3
    [t_obs, cens, T] = simulate(s, n, mu, t_cap, mu_c);
    d_true = histcounts(T, edges) / n / bin_w;
    cut_t{s} = t_obs(cens);

    ax_top(s) = nexttile(s); hold on
    bin_stairs(ax_top(s), edges, d_true, truth_col, 7, 'full sample, no cut');
    ax_dif = nexttile(3 + s); hold on
    yline(ax_dif, 0, 'Color', truth_col, 'LineWidth', 2, 'HandleVisibility', 'off');

    % omit first and thicker, since it lies under as recorded where the two agree
    for m = [2 1 3]
        [mass, covered] = bin_mass(m, t_obs, cens, edges);
        d = mass / bin_w;
        d(~covered) = NaN;
        bin_stairs(ax_top(s), edges, d, method_col(m, :), line_width(m), method_name{m});
        bin_stairs(ax_dif, edges, d - d_true, method_col(m, :), line_width(m), '');
    end
    yl = ylim(ax_dif);
    ylim(ax_dif, [-.1 .1]);

    title(ax_top(s), {scenario_name{s}, sprintf('%d of %d cut', nnz(cens), n)}, 'FontWeight', 'normal');
    xlim(ax_top(s), [0 10]);
    xlim(ax_dif, [0 10]);
    apply_generic(ax_dif, 'font_size', 16, 'tick_length', 0.025, 'xticks', [0 5 10], 'line_width', 2);
    if s == 1
        ylabel(ax_dif, 'estimate - full sample');
    end
end

% Same scale in the three top panels, with room under them for the cut freezes
y_max = max(arrayfun(@(a) a.YLim(2), ax_top));
for s = 1:3
    ylim(ax_top(s), [-0.06  1.06]);
    plot(ax_top(s), cut_t{s}, -0.035 * y_max * ones(size(cut_t{s})), '|', 'Color', [0.8 0.1 0.1], 'MarkerSize', 10, 'HandleVisibility', 'off');
    apply_generic(ax_top(s), 'font_size', 16, 'tick_length', 0.025, 'xticks', [0 5 10], 'line_width', 2);
end
ylabel(ax_top(1), 'density (1/s)');
legend(ax_top(1), 'Location', 'northeast', 'Box', 'off', 'FontSize', 13);
linkaxes(ax_top, 'y')
xlabel(tl, 'Freeze duration (s)', 'FontSize', 20);

%% The weights worked out, around the cut freeze of A
rng(1)
[t_obs, cens] = simulate(1, n, mu, t_cap, mu_c);
mass_km = bin_mass(3, t_obs, cens, edges);
ended   = histcounts(t_obs(~cens), edges);          % freezes seen to end in each bin
k       = find(edges(1:end-1) <= t_cap, 1, 'last'); % the bin that starts at the cut
rows    = max(1, k - 3):min(numel(mass_km), k + 2);
plain   = ended(rows)' / n;
tab_w   = table(edges(rows)', edges(rows + 1)', ended(rows)', plain, mass_km(rows)', mass_km(rows)' ./ plain, ...
    'VariableNames', {'from', 'to', 'ended', 'plain', 'KM', 'ratio'});
fprintf('\nA, bin masses around the freeze cut at %g s\n', t_cap);
fprintf('plain = ended / n, every freeze weighing 1/n; KM = Kaplan-Meier mass of the bin; ratio = KM / plain\n');
disp(tab_w)

%% By hand against ecdf
rng(1)
fprintf('\nlargest difference between the bin masses summed from the Kaplan-Meier product and those of the ecdf steps\n');
for s = 1:3
    [t_obs, cens] = simulate(s, n, mu, t_cap, mu_c);
    tab = km_by_hand(t_obs, cens);
    b = discretize(tab.t, edges);
    ok = ~isnan(b);
    by_hand = accumarray(b(ok), tab.dF(ok), [numel(edges) - 1, 1])';
    fprintf('  %-40s %.2e\n', scenario_name{s}, max(abs(by_hand - bin_mass(3, t_obs, cens, edges))));
end

%% Error of each estimate over many samples
bin_eval = [1 3 6];                                  % s, where the bins of the table start
j_eval   = round(bin_eval / bin_w) + 1;
p_exact  = exp(-edges(1:end-1) / mu) - exp(-edges(2:end) / mu);
hist_name = [method_name, {'no cut'}];
est_bin = nan(n_rep, 4, numel(j_eval));
err_all = nan(n_rep, 4);
fprintf('\nmean of (estimate - exact bin mass) over %d samples of %d freezes, in brackets as %% of the exact mass\n', n_rep, n);
fprintf('last column: mean |estimate - exact| over the bins that are defined, which includes the sampling noise\n');
fprintf('%-26s%s\n', 'exact mass', sprintf('%17.4f', p_exact(j_eval)));
fprintf('%-26s', '');
for j = 1:numel(bin_eval)
    fprintf('%17s', sprintf('%g-%g s', bin_eval(j), bin_eval(j) + bin_w));
end
fprintf('   %s\n', 'mean |err|');
for s = 1:3
    for r = 1:n_rep
        [t_obs, cens, T] = simulate(s, n, mu, t_cap, mu_c);
        for m = 1:4
            if m <= 3
                [mass, covered] = bin_mass(m, t_obs, cens, edges);
                mass(~covered) = NaN;
            else
                mass = histcounts(T, edges) / n;
            end
            est_bin(r, m, :) = mass(j_eval);
            err_all(r, m) = mean(abs(mass - p_exact), 'omitnan');
        end
    end
    bias = squeeze(mean(est_bin, 1, 'omitnan')) - p_exact(j_eval);
    fprintf('%s\n', scenario_name{s});
    for m = 1:4
        fprintf('  %-24s', hist_name{m});
        for j = 1:numel(j_eval)
            if isnan(bias(m, j))
                fprintf('%17s', 'n/a');
            else
                fprintf('  %+.4f (%+4.0f%%)', bias(m, j), 100 * bias(m, j) / p_exact(j_eval(j)));
            end
        end
        fprintf('   %.4f\n', mean(err_all(:, m)));
    end
end
fprintf('n/a: the bin reaches past the longest observation, which was cut, so Kaplan-Meier has no mass for it\n');

%% =====================================================================
%% ============================ HELPERS ================================
%% =====================================================================

function [t_obs, cens, T] = simulate(s, n, mu, t_cap, mu_c)
% Observed times and cut flags of one sample of n freezes. T is the true duration
T = -mu * log(rand(n, 1));
switch s
    case 1
        past = find(T > t_cap);
        cens = false(n, 1);
        cens(past(randi(numel(past)))) = true;
        t_obs = T;
        t_obs(cens) = t_cap;
    case 2
        cens = T > t_cap;
        t_obs = min(T, t_cap);
    case 3
        C = -mu_c * log(rand(n, 1));
        cens = C < T;
        t_obs = min(T, C);
end
end

% ---------------------------------------------------------------------

function [x, f, last] = estimate(method, t_obs, cens)
% Step function of one treatment of the cut freezes, as ecdf returns it. last is the
% time the curve is defined up to: later than x(end) when the longest observation was
% cut, Inf otherwise
switch method
    case 1
        [f, x] = ecdf(t_obs);
    case 2
        [f, x] = ecdf(t_obs(~cens));
    case 3
        [f, x] = ecdf(t_obs, 'Censoring', cens);
end
last = Inf;
if f(end) < 1 - 1e-9
    last = max(t_obs);
end
end

% ---------------------------------------------------------------------

function F = eval_cdf(x, f, t)
% Value of the step function (x, f) at the times t, 0 before its first point
F = zeros(size(t));
for k = 1:numel(t)
    i = find(x <= t(k), 1, 'last');
    if ~isempty(i)
        F(k) = f(i);
    end
end
end

% ---------------------------------------------------------------------

function tab = km_by_hand(t_obs, cens)
% The Kaplan-Meier product, one row per distinct observed time: n_risk freezes under
% observation just before it, d that ended there, c that were cut there, S the
% survival after it, F = 1 - S, dF the step of F
t = unique(t_obs);
n_risk = zeros(size(t));
d = zeros(size(t));
c = zeros(size(t));
for k = 1:numel(t)
    n_risk(k) = sum(t_obs >= t(k));
    d(k) = sum(t_obs == t(k) & ~cens);
    c(k) = sum(t_obs == t(k) & cens);
end
S = cumprod(1 - d ./ n_risk);
dF = [1; S(1:end-1)] .* d ./ n_risk;
tab = table(t, n_risk, d, c, S, 1 - S, dF, 'VariableNames', {'t', 'n_risk', 'd', 'c', 'S', 'F', 'dF'});
end

% ---------------------------------------------------------------------

function [mass, covered] = bin_mass(method, t_obs, cens, edges)
% Probability mass in every bin of one treatment of the cut freezes. Kaplan-Meier
% takes it from the steps of its curve. covered is false for the bins that reach past
% the longest observation when that one was cut: their Kaplan-Meier mass is not known
n = numel(t_obs);
covered = true(1, numel(edges) - 1);
switch method
    case 1
        mass = histcounts(t_obs, edges) / n;
    case 2
        mass = histcounts(t_obs(~cens), edges) / nnz(~cens);
    case 3
        [f, x] = ecdf(t_obs, 'Censoring', cens);
        step = diff(f);
        b = discretize(x(2:end), edges);
        ok = ~isnan(b);
        mass = accumarray(b(ok), step(ok), [numel(edges) - 1, 1])';
        if f(end) < 1 - 1e-9
            covered = edges(2:end) <= max(t_obs) + 1e-9;
        end
end
end

% ---------------------------------------------------------------------

function bin_stairs(ax, edges, h, col, lw, name)
% The height of every bin as a step line, a gap where it is NaN
args = {'Color', col, 'LineWidth', lw};
if isempty(name)
    args = [args, {'HandleVisibility', 'off'}];
else
    args = [args, {'DisplayName', name}];
end
stairs(ax, edges, [h h(end)], args{:});
end
