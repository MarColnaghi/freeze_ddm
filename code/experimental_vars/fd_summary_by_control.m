function summary = fd_summary_by_control(T, vars, varargin)
% FD_SUMMARY_BY_CONTROL  The freeze duration distribution at each level of each
% control variable, in numbers.
%
%   summary = fd_summary_by_control(T, {'ls', 'moving_flies', 'nloom'})
%
% T is the table of loom-evoked freezes (loom_evoked_freezes), with durations_s,
% is_censored, fly and ls. For every variable, level and loom speed it gives the
% number of freezes and of flies, the share that is censored, and the
% Kaplan-Meier median duration and probability that a freeze lasts beyond each of
% times_s. Both are conditional on the freeze having lasted the minimum duration,
% since every freeze in T did, and contact, the end of the loom cycle and the end
% of the recording count as censoring, as in ecdf(..., 'Censoring').
%
% The 95% intervals resample flies with replacement, within loom speed, and take
% the 2.5% and 97.5% of the resampled estimates. The freezes of one fly are not
% independent, so resampling bouts would give intervals that are too narrow.
%
% Variables: 'ls' (loom speed), 'moving_flies', 'nloom' (in bins of 5 looms) or any
% numeric column, which is cut into its quartiles over all the freezes, with the
% same edges as the quartiles of fd_distr_withparam_new.
%
% Options
%   times_s   [1 3 6]   durations the survival probability is taken at
%   n_boot    200       resamples of the flies
%   seed      1
%   verbose   true      print the table

opt = inputParser;
addParameter(opt, 'times_s', [1 3 6]);
addParameter(opt, 'n_boot', 200);
addParameter(opt, 'seed', 1);
addParameter(opt, 'verbose', true);
parse(opt, varargin{:});

times_s = opt.Results.times_s(:)';
n_boot = opt.Results.n_boot;
% The seed makes the intervals repeatable; the caller's random generator is put back on exit
rng_state = rng;
restore_rng = onCleanup(@() rng(rng_state));
rng(opt.Results.seed);

t = T.durations_s;
c = logical(T.is_censored);
[~, ~, fly_idx] = unique(T.fly);
n_fly = max(fly_idx);
fly_speed = accumarray(fly_idx, T.ls, [n_fly 1], @max);
speed_names = {'slow', 'fast'};

%% Cells: one per variable, level and loom speed

cells = struct('variable', {}, 'level', {}, 'speed', {}, 'rows', {});
for v = 1:numel(vars)
    [level, level_names] = control_levels(T, vars{v});
    for k = 1:numel(level_names)
        for s = [0 1]
            if strcmp(vars{v}, 'ls') && s ~= k - 1, continue, end
            rows = find(level == k & T.ls == s);
            if isempty(rows), continue, end
            cells(end + 1) = struct('variable', vars{v}, 'level', level_names{k}, 'speed', s, 'rows', rows); %#ok<AGROW>
        end
    end
end
n_cell = numel(cells);
n_stat = 1 + numel(times_s);

%% Estimates, and the same on resampled flies

point = nan(n_cell, n_stat);
for i = 1:n_cell
    r = cells(i).rows;
    point(i, :) = km_stats(t(r), c(r), ones(numel(r), 1), times_s);
end

boot = nan(n_boot, n_cell, n_stat);
for b = 1:n_boot
    counts = zeros(n_fly, 1);
    for s = [0 1]
        members = find(fly_speed == s);
        if isempty(members), continue, end % a loom speed that was not run, LC11 has only the fast one
        counts = counts + accumarray(members(randi(numel(members), numel(members), 1)), 1, [n_fly 1]);
    end
    for i = 1:n_cell
        r = cells(i).rows;
        w = counts(fly_idx(r));
        keep = w > 0;
        boot(b, i, :) = km_stats(t(r(keep)), c(r(keep)), w(keep), times_s);
    end
end

lo = nan(n_cell, n_stat);
hi = nan(n_cell, n_stat);
for i = 1:n_cell
    for j = 1:n_stat
        x = sort(boot(~isnan(boot(:, i, j)), i, j));
        if numel(x) >= 0.8 * n_boot
            lo(i, j) = x(max(1, round(0.025 * numel(x))));
            hi(i, j) = x(min(numel(x), round(0.975 * numel(x))));
        end
    end
end

%% Table

variable = string({cells.variable})';
level = string({cells.level})';
speed = string(speed_names([cells.speed] + 1))';
n = arrayfun(@(q) numel(q.rows), cells)';
n_flies = arrayfun(@(q) numel(unique(fly_idx(q.rows))), cells)';
pct_censored = arrayfun(@(q) 100 * mean(c(q.rows)), cells)';

summary = table(variable, level, speed, n, n_flies, pct_censored, ...
    point(:, 1), lo(:, 1), hi(:, 1), 'VariableNames', {'variable', 'level', 'speed', 'n', 'n_flies', 'pct_censored', 'median_s', 'median_lo', 'median_hi'});
for j = 1:numel(times_s)
    name = matlab.lang.makeValidName(sprintf('S_%gs', times_s(j)));
    summary.(name) = point(:, 1 + j);
    summary.([name '_lo']) = lo(:, 1 + j);
    summary.([name '_hi']) = hi(:, 1 + j);
end

if opt.Results.verbose
    for v = 1:numel(vars)
        m = find(summary.variable == vars{v});
        fprintf('\n%s\n', vars{v});
        fprintf('  %-6s %-5s %6s %6s %6s  %-22s', 'level', 'speed', 'n', 'flies', 'cens%', 'median s (95% CI)');
        for j = 1:numel(times_s), fprintf('  %-22s', sprintf('P(>%g s) (95%% CI)', times_s(j))); end
        fprintf('\n');
        for i = m'
            fprintf('  %-6s %-5s %6d %6d %5.0f%%  %-22s', summary.level(i), summary.speed(i), n(i), n_flies(i), pct_censored(i), fmt_ci(point(i, 1), lo(i, 1), hi(i, 1), '%.2f'));
            for j = 1:numel(times_s), fprintf('  %-22s', fmt_ci(point(i, 1 + j), lo(i, 1 + j), hi(i, 1 + j), '%.2f')); end
            fprintf('\n');
        end
    end
end
end

%% Local functions

function out = km_stats(t, c, w, times_s)
% Kaplan-Meier median and survival at times_s, from durations t with censoring c and
% frequencies w. NaN where the estimate is not defined: the curve never reaches one
% half, or the follow-up ends before the time asked for
out = nan(1, 1 + numel(times_s));
if isempty(t) || all(c), return, end
[F, x] = ecdf(t, 'Censoring', c, 'Frequency', w);
k = find(F >= 0.5, 1);
if ~isempty(k), out(1) = x(k); end
for j = 1:numel(times_s)
    if times_s(j) > max(t), continue, end
    i = find(x <= times_s(j), 1, 'last');
    if isempty(i), out(1 + j) = 1; else, out(1 + j) = 1 - F(i); end
end
end

function [level, names] = control_levels(T, var)
switch var
    case 'ls'
        level = T.ls + 1;
        names = {'slow', 'fast'};
    case 'late_onset'
        level = T.late_onset + 1;
        names = {'early', 'late'};
    case 'moving_flies'
        level = T.moving_flies + 1;
        names = arrayfun(@num2str, 0:4, 'UniformOutput', false);
    case {'nloom', 'ln'}
        level = discretize(T.nloom, 0.5:5:20.5);
        names = {sprintf('%d-5', max(1, min(T.nloom))), '6-10', '11-15', '16-20'};
    otherwise
        % The same quartile edges as the figure (fd_distr_withparam_new), so that
        % Q1-Q4 are the same freezes in both
        edges = unique(quantile(T.(var), [0 0.25 0.5 0.75 1]));
        level = discretize(T.(var), edges);
        names = arrayfun(@(q) sprintf('Q%d', q), 1:numel(edges) - 1, 'UniformOutput', false);
end
end

function s = fmt_ci(x, lo, hi, f)
if isnan(x)
    s = 'not reached';
elseif isnan(lo)
    s = sprintf(f, x);
else
    s = sprintf([f ' (' f '-' f ')'], x, lo, hi);
end
end
