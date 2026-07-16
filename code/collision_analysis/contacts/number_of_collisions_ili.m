clearvars; close all; clc;

% Colors + paths
col   = cmapper([], 5);
paths = path_generator('folder', 'collision_analysis/contacts');

% Shared data loading (see load_contact_timeseries)
ls   = 25;                       % loom speed
data = load_contact_timeseries(ls);

md_mat         = data.md;
n_moving_flies = data.n_moving_flies;
n_flies        = data.n_flies;

% Loom onset frames (+ end-of-recording sentinel for windowed slicing)
loom_times = [data.loom_times, 36000 * ones(n_flies, 1)];
nIntervals = data.n_looms;

% Min-distance slices in a window after each loom onset
window    = 660;
md_slices = arrayfun(@(fly, k) ...
    md_mat(fly, loom_times(fly,k) : min(loom_times(fly,k)+window, size(md_mat,2))), ...
    repmat((1:n_flies)', 1, nIntervals), ...
    repmat(1:nIntervals,  n_flies,    1), ...
    'UniformOutput', false);

% Contact detection
threshold_contact = 70;
contact = cellfun(@(x) any(x < threshold_contact), md_slices);

% --- Contact events per fly (sanity diagnostic) ---
fh = figure('Color', 'w', 'Position', [100 100 560 360]);
plot(sum(contact, 2), 'Color', col.var.mindist, 'LineWidth', 2)
apply_generic(gca, 'font_size', 20)
xlabel('Fly index')
ylabel('Contacts for each Loom')
exporter(fh, paths, 'contacts_per_fly.pdf')

% --- First-contact latency ---
first_contact_frame          = NaN(size(md_slices));
first_contact_frame(contact) = cellfun( ...
    @(x) find(x < threshold_contact, 1, 'first'), ...
    md_slices(contact));
latency_seconds = first_contact_frame / 60;

% Control: same procedure on random baseline windows (frames 1-18000).
% Each fly contributes nIntervals random start times — same N as the loom condition.
baseline_end    = min(18000, size(md_mat,2));
baseline_margin = baseline_end - window;          % latest valid start frame
rng(42);                                          % reproducibility
first_contact_baseline = NaN(n_flies, nIntervals);
for f = 1:n_flies
    t_starts = sort(randperm(baseline_margin, nIntervals));
    for k = 1:nIntervals
        seg = md_mat(f, t_starts(k) : t_starts(k)+window);
        hit = find(seg < threshold_contact, 1, 'first');
        if ~isempty(hit)
            first_contact_baseline(f, k) = hit;
        end
    end
end
latency_baseline_s = first_contact_baseline / 60;

% Latency distribution: keep only contacts that FORM after the loom. Pairs already
% in contact at onset (first frame below threshold -> latency == 1 frame) are
% pre-existing proximity, not loom-evoked approaches, so exclude them from the
% histogram. This affects ONLY the probabilities plotted below; the printed contact
% rates and every other panel still use all contacts.
latency_loom_formed = latency_seconds;
latency_loom_formed(first_contact_frame == 1) = NaN;
latency_base_formed = latency_baseline_s;
latency_base_formed(first_contact_baseline == 1) = NaN;

edges = (0:10/60:10.5) - 1/120;
fh = figure('Color', 'w', 'Position', [100 100 6300 360]);
hold on
histogram(latency_loom_formed(:), edges, 'Normalization', 'probability', ...
    'FaceColor', hex2rgb(col.period.loom), 'FaceAlpha', 0.6, 'EdgeColor', 'none')
histogram(latency_base_formed(:), edges, 'Normalization', 'probability', ...
    'FaceColor', hex2rgb(col.period.bsl),  'FaceAlpha', 0.6, 'EdgeColor', 'none')
apply_generic(gca, 'font_size', 20)
xlabel('Latency to contact (s)')
ylabel('Probability')
legend({'Loom-triggered', 'Baseline (random windows)'}, 'box', 'off', 'FontSize', 16)
% Report the fraction of pairs excluded above (already in contact at onset)
frac_touching_loom = 100 * mean(first_contact_frame(:)    == 1);
frac_touching_base = 100 * mean(first_contact_baseline(:) == 1);
title(sprintf('Already touching at onset — loom %.0f%%, baseline %.0f%%', ...
    frac_touching_loom, frac_touching_base), 'FontWeight', 'normal', 'FontSize', 14)
exporter(fh, paths, 'contact_latency_loom_vs_baseline.pdf')

fprintf('Baseline contact rate:     %.1f%%\n', 100*mean(~isnan(latency_baseline_s(:))))
fprintf('Loom contact rate:         %.1f%%\n', 100*mean(~isnan(latency_seconds(:))))

% --- Composition of pairs: already touching / contact forms / no contact ---
% Puts the 'already at onset' fraction (previously only in the histogram title)
% in proportion against the formed-contact and no-contact fractions. Loom and
% baseline share the same N = n_flies * nIntervals, so the bars are comparable.
%   first_contact_* : NaN  -> no contact in window
%                     == 1  -> already below threshold at onset (excluded above)
%                     >  1  -> contact forms during the window
categorize = @(m) [mean(m(:) == 1), mean(m(:) > 1), mean(isnan(m(:)))];
comp   = [categorize(first_contact_frame); categorize(first_contact_baseline)];
counts = [sum(first_contact_frame(:)    == 1), sum(first_contact_frame(:)    > 1), sum(isnan(first_contact_frame(:))); ...
          sum(first_contact_baseline(:) == 1), sum(first_contact_baseline(:) > 1), sum(isnan(first_contact_baseline(:)))];

cat_col = [0.85 0.33 0.10;    % already touching at onset
           0.20 0.55 0.75;    % contact forms in window
           0.85 0.85 0.85];   % no contact

fh = figure('Color', 'w', 'Position', [100 100 560 420]);
b = bar(100*comp, 'stacked', 'EdgeColor', 'none');
for s = 1:3
    b(s).FaceColor = cat_col(s, :);
end
yctr = 100*cumsum(comp, 2) - 100*comp/2;          % centre of each stacked segment
for g = 1:2
    for s = 1:3
        if comp(g, s) > 0.01
            text(g, yctr(g, s), sprintf('%.0f%%\n(n=%d)', 100*comp(g, s), counts(g, s)), ...
                'HorizontalAlignment', 'center', 'FontSize', 12)
        end
    end
end
apply_generic(gca, 'font_size', 20, 'xticks', [1 2])
xticklabels({'Loom', 'Baseline'})
ylabel('% of (fly, loom) pairs')
ylim([0 100])
legend({'Already touching at onset', 'Contact forms in window', 'No contact'}, ...
    'box', 'off', 'FontSize', 12, 'Location', 'eastoutside')
exporter(fh, paths, 'contact_composition_loom_vs_baseline.pdf')

% --- Distance at and just before loom onset (diagnostics) ---
distance_at_loom = arrayfun(@(fly, k) ...
    md_mat(fly, loom_times(fly,k)), ...
    repmat((1:n_flies)', 1, nIntervals), ...
    repmat(1:nIntervals,  n_flies,    1));
distance_before_loom = arrayfun(@(fly, k) ...
    md_mat(fly, max(loom_times(fly,k)-1, 1)), ...   % guard against frame-1 loom onset
    repmat((1:n_flies)', 1, nIntervals), ...
    repmat(1:nIntervals,  n_flies,    1));
n_already  = sum(distance_at_loom(:)    < threshold_contact);
n_latency0 = sum(first_contact_frame(:) == 1, 'omitnan');
fprintf('Already contacting at loom onset:    %d  (%.1f%%)\n', n_already, 100*mean(distance_at_loom(:) < threshold_contact));
fprintf('Already contacting one frame before: %.1f%%\n',                  100*mean(distance_before_loom(:) < threshold_contact));
fprintf('Contacts with latency = 1 frame:     %d\n', n_latency0);

% --- Loom-triggered pooled mean min-distance (overlaid on the per-condition plot below) ---
% Pre-onset window is critical — if distance is already decreasing before t=0,
% the dip is NOT caused by the loom (flies were already approaching).
pre  = 180;   % 3 s before loom onset (60 fps)
post = 630;   % 10.5 s after loom onset (60 fps)
t_axis   = -pre:post;
traj_all = NaN(n_flies * nIntervals, pre+post+1);
pair = 0;
for f = 1:n_flies
    for k = 1:nIntervals
        pair    = pair + 1;
        t0      = loom_times(f, k);
        t_start = max(t0 - pre,  1);
        t_end   = min(t0 + post, size(md_mat, 2));
        traj_all(pair, (t_start - t0 + pre + 1):(t_end - t0 + pre + 1)) = md_mat(f, t_start:t_end);
    end
end
mu = mean(traj_all, 1, 'omitnan');   % pooled mean across all (fly, loom) pairs

% --- Contact rate by condition (number of moving conspecifics) ---
conditions   = unique(n_moving_flies)';
contact_rate = mean(contact, 2);                  % fraction of looms with a contact, per fly
group_means  = NaN(1, numel(conditions));
group_sems   = NaN(1, numel(conditions));
for c = 1:numel(conditions)
    rates          = contact_rate(n_moving_flies == conditions(c));
    group_means(c) = mean(rates);
    group_sems(c)  = std(rates) / sqrt(numel(rates));
end

fh = figure('Color', 'w', 'Position', [100 100 520 420]);
hold on
for c = 1:numel(conditions)
    idx = n_moving_flies == conditions(c);
    scatter(conditions(c) + 0.05*randn(sum(idx),1), contact_rate(idx), 20, ...
        col.n_mov_flies(c,:), 'filled', 'MarkerFaceAlpha', 0.5)
end
errorbar(conditions, group_means, group_sems, 'ko-', 'LineWidth', 1.5, 'MarkerFaceColor', 'k')
apply_generic(gca, 'font_size', 20, 'xticks', conditions)
xlabel('Number of Moving Flies')
ylabel({'Fraction of Looms', 'with Contact'})
exporter(fh, paths, 'contact_rate_by_condition.pdf')

% --- Loom-triggered average min-distance per condition ---
fh = figure('Color', 'w', 'Position', [100 100 560 360]);
hold on
for c = 1:numel(conditions)
    flies_c = find(n_moving_flies == conditions(c));
    pairs_c = NaN(numel(flies_c) * nIntervals, pre+post+1);
    row = 0;
    for fi = 1:numel(flies_c)
        f = flies_c(fi);
        for k = 1:nIntervals
            row     = row + 1;
            t0      = loom_times(f, k);
            t_start = max(t0 - pre,  1);
            t_end   = min(t0 + post, size(md_mat, 2));
            pairs_c(row, (t_start - t0 + pre + 1):(t_end - t0 + pre + 1)) = md_mat(f, t_start:t_end);
        end
    end
    plot(t_axis, mean(pairs_c, 1, 'omitnan'), 'Color', col.n_mov_flies(c,:), 'LineWidth', 2)
end
plot(t_axis, mu, 'k', 'LineWidth', 2.5)                       % pooled mean across all conditions
xline(0, '--', 'Color', 'k', 'LineWidth', 1, 'HandleVisibility', 'off')
yline(threshold_contact, ':', 'Color', 'k', 'LineWidth', 1, 'HandleVisibility', 'off')
apply_generic(gca, 'font_size', 20)
xlabel('Frames relative to loom onset')
ylabel('Mean min inter-fly distance')
legend([arrayfun(@(c) sprintf('%d moving', c), conditions, 'UniformOutput', false), {'all'}], ...
    'box', 'off', 'FontSize', 14)
exporter(fh, paths, 'loom_triggered_distance_by_condition.pdf')

% --- Median min-distance over the full timeline by condition ---
median_loom_ts = median(loom_times(:, 1:nIntervals), 1);
fh = figure('Color', 'w', 'Position', [100 100 700 360]);
hold on
for c = 1:numel(conditions)
    plot(median(md_mat(n_moving_flies == conditions(c), :), 1), 'Color', col.n_mov_flies(c,:), 'LineWidth', 1.5)
end
xline(median_loom_ts, 'Color', [0.7 0.7 0.7])
apply_generic(gca, 'font_size', 20)
xlabel('Frame')
ylabel('Median min inter-fly distance')
legend(arrayfun(@(c) sprintf('%d moving', c), conditions, 'UniformOutput', false), 'box', 'off', 'FontSize', 14)
exporter(fh, paths, 'median_distance_timeline_by_condition.pdf')
