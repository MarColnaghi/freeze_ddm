% Split-half check of the loom-evoked definition. The windows and the minimum
% duration were chosen by looking at every fly and then scored on the same flies,
% which flatters them. Here the flies are split at random into two halves, the
% definition is learned on one half and scored on the other, and that is repeated
% over many splits, to see how far the edges move from split to split and how
% much purity and retention drop on flies the choice never saw.
%
% What is learned on a half of the flies, for each loom speed:
%   core        the run of latencies around the main peak where at least
%               core_fraction of the onsets are in excess of the baseline mock
%               looms (3-frame bins)
%   inclusive   the narrowest interval holding inclusive_share of the excess
%   minimum duration   the shortest one in min_durations for which the core
%               window has a purity of at least purity_target
% Purity is 1 - baseline / loom onsets inside a window. Kept is the excess inside
% it over the whole excess (latency -10 to 150, at the minimum duration of
% loom_evoked_bouts) of the same flies, so it counts the share of all evoked onsets
% a definition holds on to.
%
% The split is by fly, since the bouts of one fly are not independent and each fly
% is one recording, and it is balanced over loom speed and number of moving flies.
% A half holds half the data, so the edges move more here than they would for the
% full sample.

[bouts, thresholds] = loom_evoked_data;

n_splits = 200;                         % random splits of the flies into two halves
purity_target = 0.90;
core_fraction = 0.75;
inclusive_share = 0.90;
first_min = thresholds.min_dur;         % frames, the minimum of loom_evoked_bouts
min_durations = [first_min:30 35 40 50 60];
fixed_durations = [first_min 15 20 30]; % minimum durations the fixed windows are scored at
latency = -10:150;                      % frames after the loom the excess is looked at
speeds = [25 50];
rng(1);

% Windows that are not learned, scored on the same halves. Rows are the loom
% speeds, columns the windows
fixed_names = {'under test', '[5 55]'};
fixed_windows = {thresholds.le_window_sl, [5 55]; thresholds.le_window_fl, [5 55]};

bouts_loom = loom_evoked_bouts(bouts);
bouts_base = loom_evoked_bouts(bouts, 'period', 'bsl');

%% Onset counts per fly, latency and minimum duration

flies = unique([bouts_loom.fly; bouts_base.fly]);
n_flies = numel(flies);
[~, loom_fly] = ismember(bouts_loom.fly, flies);
[~, base_fly] = ismember(bouts_base.fly, flies);

fly_speed = nan(n_flies, 1);
fly_moving = nan(n_flies, 1);
fly_speed(loom_fly) = bouts_loom.sloom;
fly_moving(loom_fly) = bouts_loom.moving_flies;
fly_speed(base_fly) = bouts_base.sloom;
fly_moving(base_fly) = bouts_base.moving_flies;
assert(~any(isnan(fly_speed)) && all(fly_speed(loom_fly) == bouts_loom.sloom) && all(fly_speed(base_fly) == bouts_base.sloom), ...
    'split_half_check:flySpeed', 'A fly has bouts at two loom speeds.')

count_loom = onset_counts(bouts_loom, loom_fly, n_flies, latency, min_durations);
count_base = onset_counts(bouts_base, base_fly, n_flies, latency, min_durations);

fixed_index = arrayfun(@(d) find(min_durations == d), fixed_durations);
strata = fly_speed * 10 + fly_moving;
stratum_ids = unique(strata);

%% Splits

n_pairs = 2 * n_splits;          % each split is learned on A and scored on B, and the other way round
n_fixed = numel(fixed_names);
n_defs = 4;                      % core and inclusive, at the first minimum duration and at the learned one

for si = 1:numel(speeds)
    res(si).core = nan(n_pairs, 2);
    res(si).inclusive = nan(n_pairs, 2);
    res(si).min_duration = nan(n_pairs, 1);
    res(si).purity_in = nan(n_pairs, n_defs);
    res(si).purity_out = nan(n_pairs, n_defs);
    res(si).kept_in = nan(n_pairs, n_defs);
    res(si).kept_out = nan(n_pairs, n_defs);
    res(si).fixed_purity = nan(n_pairs, n_fixed, numel(fixed_durations));
    res(si).fixed_kept = nan(n_pairs, n_fixed, numel(fixed_durations));
end

for idx_split = 1:n_splits
    in_a = false(n_flies, 1);
    for idx_stratum = 1:numel(stratum_ids)
        members = find(strata == stratum_ids(idx_stratum));
        members = members(randperm(numel(members)));
        in_a(members(1:ceil(numel(members) / 2))) = true;
    end
    halves = {find(in_a), find(~in_a)};

    for si = 1:numel(speeds)
        of_speed = fly_speed == speeds(si);

        for h = 1:2
            train = halves{h}(of_speed(halves{h}));
            test = halves{3 - h}(of_speed(halves{3 - h}));
            pair = (idx_split - 1) * 2 + h;

            [l_train, b_train] = half_counts(count_loom, count_base, train);
            [l_test, b_test] = half_counts(count_loom, count_base, test);

            % Learn on the training half
            [core, inclusive] = learn_windows(l_train(:, 1), b_train(:, 1), latency, core_fraction, inclusive_share);
            k_min = learn_min_duration(l_train, b_train, latency, core, purity_target);
            res(si).core(pair, :) = core;
            res(si).inclusive(pair, :) = inclusive;
            res(si).min_duration(pair) = min_durations(k_min);

            % Score it on the training half and on the other half
            windows = {core, inclusive, core, inclusive};
            columns = [1 1 k_min k_min];
            for d = 1:n_defs
                [res(si).purity_in(pair, d), res(si).kept_in(pair, d)] = score_window(l_train(:, columns(d)), b_train(:, columns(d)), latency, windows{d}, sum(l_train(:, 1) - b_train(:, 1)));
                [res(si).purity_out(pair, d), res(si).kept_out(pair, d)] = score_window(l_test(:, columns(d)), b_test(:, columns(d)), latency, windows{d}, sum(l_test(:, 1) - b_test(:, 1)));
            end

            % Windows that are not learned, on this half
            for k = 1:n_fixed
                for t = 1:numel(fixed_durations)
                    [res(si).fixed_purity(pair, k, t), res(si).fixed_kept(pair, k, t)] = score_window(l_train(:, fixed_index(t)), b_train(:, fixed_index(t)), latency, fixed_windows{si, k}, sum(l_train(:, 1) - b_train(:, 1)));
                end
            end
        end
    end
end

%% Report

def_names = {'core, first min', 'inclusive, first min', 'core, learned min', 'inclusive, learned min'};

for si = 1:numel(speeds)
    r = res(si);
    all_flies = find(fly_speed == speeds(si));
    [l_all, b_all] = half_counts(count_loom, count_base, all_flies);
    [core_all, inclusive_all] = learn_windows(l_all(:, 1), b_all(:, 1), latency, core_fraction, inclusive_share);
    min_all = min_durations(learn_min_duration(l_all, b_all, latency, core_all, purity_target));

    fprintf('\n%d cm/s: %d flies, %d per half, %d splits\n', speeds(si), numel(all_flies), floor(numel(all_flies) / 2), n_splits);
    fprintf('  all flies:   core %s | inclusive %s | minimum duration %d\n', mat2str(core_all), mat2str(inclusive_all), min_all);
    fprintf('  learned on half the flies, median (2.5%% to 97.5%%) over %d halves:\n', n_pairs);
    fprintf('    core start %g (%g to %g) | core end %g (%g to %g)\n', pctl(r.core(:, 1), 0.5), pctl(r.core(:, 1), 0.025), pctl(r.core(:, 1), 0.975), pctl(r.core(:, 2), 0.5), pctl(r.core(:, 2), 0.025), pctl(r.core(:, 2), 0.975));
    fprintf('    inclusive start %g (%g to %g) | inclusive end %g (%g to %g)\n', pctl(r.inclusive(:, 1), 0.5), pctl(r.inclusive(:, 1), 0.025), pctl(r.inclusive(:, 1), 0.975), pctl(r.inclusive(:, 2), 0.5), pctl(r.inclusive(:, 2), 0.025), pctl(r.inclusive(:, 2), 0.975));
    fprintf('    minimum duration for a core purity of %.2f: %g (%g to %g)\n', purity_target, pctl(r.min_duration, 0.5), pctl(r.min_duration, 0.025), pctl(r.min_duration, 0.975));

    fprintf('  learned on one half, scored on the same half and on the other (mean over %d, in -> out):\n', n_pairs);
    fprintf('    %-24s %-18s %s\n', 'definition', 'purity', 'kept');
    for d = 1:n_defs
        fprintf('    %-24s %.3f -> %.3f      %3.0f%% -> %3.0f%%\n', def_names{d}, mean(r.purity_in(:, d)), mean(r.purity_out(:, d)), 100 * mean(r.kept_in(:, d)), 100 * mean(r.kept_out(:, d)));
    end

    fprintf('  windows that are not learned, over %d halves (mean, 2.5%% to 97.5%%):\n', n_pairs);
    fprintf('    %-12s %-10s %-4s %-26s %s\n', 'window', 'latency', 'min', 'purity', 'kept');
    for k = 1:n_fixed
        for t = 1:numel(fixed_durations)
            p = r.fixed_purity(:, k, t);
            kp = r.fixed_kept(:, k, t);
            fprintf('    %-12s %-10s %-4d %.3f (%.3f to %.3f)      %3.0f%% (%3.0f%% to %3.0f%%)\n', fixed_names{k}, mat2str(fixed_windows{si, k}), fixed_durations(t), mean(p), pctl(p, 0.025), pctl(p, 0.975), 100 * mean(kp), 100 * pctl(kp, 0.025), 100 * pctl(kp, 0.975));
        end
    end
end

%% Figure

col = cmapper([], 5);
fh = figure('color', 'w', 'Position', [100, 100, 1500, 450]);
tl = tiledlayout(fh, 1, 3, 'TileSpacing', 'compact');
edge_names = {'core start', 'core end', 'inclusive start', 'inclusive end'};

for si = 1:numel(speeds)
    ax(si) = nexttile;
    hold on
    all_flies = find(fly_speed == speeds(si));
    [l_all, b_all] = half_counts(count_loom, count_base, all_flies);
    [core_all, inclusive_all] = learn_windows(l_all(:, 1), b_all(:, 1), latency, core_fraction, inclusive_share);
    edges_all = [core_all, inclusive_all];
    edges_split = [res(si).core, res(si).inclusive];
    c = col.vars.sloom(3 * si, :);

    for e = 1:4
        plot([pctl(edges_split(:, e), 0.025), pctl(edges_split(:, e), 0.975)], [e e], '-', 'Color', c, 'LineWidth', 6);
        plot(pctl(edges_split(:, e), 0.5), e, 'o', 'MarkerFaceColor', 'w', 'MarkerEdgeColor', 'k', 'MarkerSize', 9, 'LineWidth', 1.5);
        plot(edges_all(e), e, 'kx', 'MarkerSize', 10, 'LineWidth', 2);
    end

    set(ax(si), 'YDir', 'reverse', 'YTick', 1:4, 'YTickLabel', edge_names, 'YLim', [0.5 4.5]);
    xlabel('Latency after the loom (frames)')
    apply_generic(ax(si), 'font_size', 16);
    title(sprintf('%d cm/s', speeds(si)), 'FontSize', 18)
end
linkaxes(ax(1:2), 'x')

ax(3) = nexttile;
hold on
for si = 1:numel(speeds)
    histogram(res(si).min_duration, 'BinMethod', 'integers', 'FaceColor', col.vars.sloom(3 * si, :), 'FaceAlpha', 0.7, 'EdgeColor', 'none');
end
xlabel(sprintf('Minimum duration reaching a core purity of %.2f (frames)', purity_target))
ylabel('Halves')
apply_generic(ax(3), 'font_size', 16);
legend({'25 cm/s', '50 cm/s'}, 'Box', 'off', 'Location', 'northeast', 'FontSize', 16)

%% Local functions

function counts = onset_counts(bouts_table, fly_index, n_flies, latency, min_durations)
% Onsets per fly, latency and minimum duration, the durations being the time each
% bout was observed for
lat_index = bouts_table.onsets_loomwin - latency(1) + 1;
valid = lat_index >= 1 & lat_index <= numel(latency);
counts = zeros(n_flies, numel(latency), numel(min_durations));
for k = 1:numel(min_durations)
    keep = valid & bouts_table.durations >= min_durations(k);
    counts(:, :, k) = accumarray([fly_index(keep), lat_index(keep)], 1, [n_flies, numel(latency)]);
end
end

function [l, b] = half_counts(count_loom, count_base, flies)
% Onsets per latency and minimum duration summed over some flies
l = reshape(sum(count_loom(flies, :, :), 1), size(count_loom, 2), size(count_loom, 3));
b = reshape(sum(count_base(flies, :, :), 1), size(count_base, 2), size(count_base, 3));
end

function [core, inclusive] = learn_windows(l, b, latency, core_fraction, inclusive_share)
% l and b are the loom and baseline onsets per latency
e = l - b;

% Core: 3-frame bins, the run around the peak where enough of the onsets are in excess
n3 = floor(numel(latency) / 3);
l3 = sum(reshape(l(1:3 * n3), 3, n3), 1);
b3 = sum(reshape(b(1:3 * n3), 3, n3), 1);
excess3 = l3 - b3;
[~, peak] = max(excess3);
enough = excess3 ./ max(l3, 1) >= core_fraction;
first = peak;
last = peak;
if enough(peak)
    while first > 1 && enough(first - 1), first = first - 1; end
    while last < n3 && enough(last + 1), last = last + 1; end
end
core = [latency(3 * (first - 1) + 1), latency(3 * last)];

% Inclusive: the narrowest interval holding inclusive_share of the net excess
total = sum(e);
best_width = inf;
best_start = 1;
for i = 1:numel(e)
    j = find(cumsum(e(i:end)) >= inclusive_share * total, 1);
    if ~isempty(j) && j < best_width
        best_width = j;
        best_start = i;
    end
end
inclusive = [latency(best_start), latency(best_start + best_width - 1)];
end

function k = learn_min_duration(l, b, latency, window, purity_target)
% Index of the shortest minimum duration whose purity in the window reaches the target
i = find(latency == window(1));
j = find(latency == window(2));
purity = 1 - sum(b(i:j, :), 1) ./ sum(l(i:j, :), 1);
k = find(purity >= purity_target, 1);
if isempty(k), k = numel(purity); end
end

function [purity, kept] = score_window(l, b, latency, window, total)
i = find(latency == window(1));
j = find(latency == window(2));
n_loom = sum(l(i:j));
n_base = sum(b(i:j));
purity = 1 - n_base / n_loom;
kept = (n_loom - n_base) / total;
end

function v = pctl(x, p)
x = sort(x(~isnan(x)));
v = x(min(numel(x), max(1, round(p * numel(x)))));
end
