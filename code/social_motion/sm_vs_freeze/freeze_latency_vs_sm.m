clearvars; close all; clc;

thresholds = define_thresholds;
paths      = path_generator('folder', 'social_motion');
bouts      = importdata(fullfile(paths.dataset, 'bouts.mat'));
bouts      = bouts_formatting(bouts, thresholds);
sm_cache   = importdata(fullfile(paths.cache_path, 'motion_cache.mat'));

% Loom-triggered freeze bouts only
% bouts.le flags onsets within the defined loom-event window (set in bouts_formatting)
sel   = bouts.le == 1 & bouts.type == 1;
bl    = bouts(sel, :);

latency_s  = bl.onsets_loomaligned / 60;   % latency from loom onset to freeze onset (s)
conditions = unique(bl.moving_flies)';      % e.g. [0 1 2 3 4]

% Social motion during the accumulation window: loom onset → freeze onset
% This is the evidence the fly accumulates before deciding to freeze.
% Not precomputed in the table, so we extract it from sm_cache here.
sm_accum = NaN(height(bl), 1);
for i = 1:height(bl)
    t_loom   = bl.loom_ts(i);
    t_freeze = bl.onsets(i);
    if t_loom > 0 && t_freeze > t_loom
        sm_fly       = sm_cache(bl.fly(i));
        sm_accum(i)  = mean(sm_fly(t_loom : t_freeze - 1));
    end
end

% ── Figure 1: CDF of latency by condition ────────────────────────────────
cmap = lines(length(conditions));
figure; hold on
for c = 1:length(conditions)
    lat_c = sort(latency_s(bl.moving_flies == conditions(c)));
    lat_c = lat_c(~isnan(lat_c));
    plot(lat_c, (1:length(lat_c)) / length(lat_c), 'Color', cmap(c,:), 'LineWidth', 1.5)
end
xlabel('Latency to freeze (s)')
ylabel('Cumulative probability')
legend(arrayfun(@(c) sprintf('%d moving', c), conditions, 'UniformOutput', false), 'Location', 'southeast')
title('Freeze latency by condition')

% ── Figure 2: Median latency ± SEM by condition ──────────────────────────
med_lat = NaN(1, length(conditions));
sem_lat = NaN(1, length(conditions));
figure; hold on
for c = 1:length(conditions)
    idx   = bl.moving_flies == conditions(c);
    lat_c = latency_s(idx);
    lat_c = lat_c(~isnan(lat_c));
    med_lat(c) = median(lat_c);
    sem_lat(c) = std(lat_c) / sqrt(length(lat_c));
    scatter(conditions(c) * ones(sum(idx), 1) + 0.08*randn(sum(idx), 1), ...
            latency_s(idx), 8, [0.8 0.8 0.8], 'filled')
end
errorbar(conditions, med_lat, sem_lat, 'ko-', 'LineWidth', 1.5, 'MarkerFaceColor', 'k')
xticks(conditions)
xlabel('Number of moving conspecifics')
ylabel('Median latency to freeze (s)')

% ── Figure 3: Social motion (loom→freeze) vs latency ────────────────────
% Core DDM prediction: more accumulation signal → shorter latency
valid = ~isnan(sm_accum) & ~isnan(latency_s);
figure; hold on
scatter(sm_accum(valid), latency_s(valid), 10, [0.7 0.7 0.7], 'filled')

% Overlay per-condition means to see the group-level trend
for c = 1:length(conditions)
    idx_c = bl.moving_flies == conditions(c) & valid;
    plot(mean(sm_accum(idx_c)), mean(latency_s(idx_c)), 'o', ...
        'MarkerFaceColor', cmap(c,:), 'MarkerEdgeColor', 'k', 'MarkerSize', 10)
end
xlabel('Mean social motion (loom onset \rightarrow freeze onset)')
ylabel('Latency to freeze (s)')
legend(['Individual bouts', arrayfun(@(c) sprintf('%d moving', c), conditions, 'UniformOutput', false)], ...
    'Location', 'northeast')

[r, p] = corr(sm_accum(valid), latency_s(valid), 'Type', 'Spearman');
title(sprintf('Spearman r = %.2f,  p = %.3f  (n = %d bouts)', r, p, sum(valid)))

fprintf('\nMedian latency by condition:\n')
for c = 1:length(conditions)
    idx = bl.moving_flies == conditions(c);
    fprintf('  %d moving: %.2f s  (n = %d bouts)\n', conditions(c), median(latency_s(idx), 'omitnan'), sum(idx))
end
fprintf('\nCorrelation sm_accum vs latency: Spearman r = %.3f, p = %.4f\n', r, p)
