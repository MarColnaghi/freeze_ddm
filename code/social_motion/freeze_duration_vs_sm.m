clearvars; close all; clc;

thresholds = define_thresholds;
paths      = path_generator('folder', 'social_motion');
bouts      = importdata(fullfile(paths.dataset, 'bouts.mat'));
bouts      = bouts_formatting(bouts, thresholds);
sm_cache   = importdata(fullfile(paths.cache_path, 'motion_cache.mat'));

% Loom-triggered freeze bouts only, both loom speeds
sel  = bouts.le == 1 & bouts.type == 1;
bl   = bouts(sel, :);

% Right-censoring: truncate duration and avg_sm at the first collision frame.
% impose_contact_threshold adds ending_time (frames) and is_censored to bl.
contact_threshold = 70;   % px — from pc_as_function_of_distance empirical estimate
bl = impose_contact_threshold(bl, 'threshold', contact_threshold, 'type', 'onlyfreeze');

% Recompute avg_sm using ending_time as the window boundary.
% For censored bouts this truncates avg_sm at the collision frame;
% for uncensored bouts it matches the original computation.
avg_sm_corrected = NaN(height(bl), 1);
for i = 1:height(bl)
    t_on  = bl.onsets(i);
    t_end = t_on + bl.ending_time(i) - 1;
    sm    = sm_cache(bl.fly(i));
    t_end = min(t_end, numel(sm));
    if t_end >= t_on
        avg_sm_corrected(i) = mean(sm(t_on:t_end));
    end
end

duration_s  = bl.ending_time / 60;   % censored duration
conditions  = unique(bl.moving_flies)';
speeds      = unique(bl.sloom)';
cmap        = lines(length(conditions));

% ── Figure 1: Duration distribution by condition (CDF) ───────────────────
for idx_sloom = 1:length(speeds)
    figure; hold on
    mask_s = bl.sloom == speeds(idx_sloom);
    for c = 1:length(conditions)
        dur_c = sort(duration_s(mask_s & bl.moving_flies == conditions(c)));
        dur_c = dur_c(~isnan(dur_c));
        plot(dur_c, (1:length(dur_c)) / length(dur_c), 'Color', cmap(c,:), 'LineWidth', 1.5)
    end
    xlabel('Freeze duration (s)')
    ylabel('Cumulative probability')
    legend(arrayfun(@(c) sprintf('%d moving', c), conditions, 'UniformOutput', false), 'Location', 'southeast')
    title(sprintf('Freeze duration — loom speed %d cm/s', speeds(idx_sloom)))
end

% ── Figure 2: Median duration ± SEM by condition ─────────────────────────
figure; hold on
for idx_sloom = 1:length(speeds)
    mask_s    = bl.sloom == speeds(idx_sloom);
    med_dur   = NaN(1, length(conditions));
    sem_dur   = NaN(1, length(conditions));
    for c = 1:length(conditions)
        idx     = mask_s & bl.moving_flies == conditions(c);
        dur_c   = duration_s(idx);
        dur_c   = dur_c(~isnan(dur_c));
        med_dur(c) = median(dur_c);
        sem_dur(c) = std(dur_c) / sqrt(length(dur_c));
    end
    errorbar(conditions + 0.1*(idx_sloom-1), med_dur, sem_dur, 'o-', ...
        'LineWidth', 1.5, 'MarkerFaceColor', [idx_sloom==2, 0, idx_sloom==1])
end
xticks(conditions)
xlabel('Number of moving conspecifics')
ylabel('Median freeze duration (s)')
legend(arrayfun(@(s) sprintf('%d cm/s', s), speeds, 'UniformOutput', false))

% ── Figure 3: avg_sm (corrected) vs censored duration ────────────────────
% Core test: does more social motion during the freeze predict shorter duration?
% Expect negative correlation if more motion = faster termination.
% Censored bouts (collision before freeze ended) are shown as open markers —
% both their duration and avg_sm are truncated at the same collision frame,
% so the relationship is internally consistent but they are right-censored.
for idx_sloom = 1:length(speeds)
    mask_s      = bl.sloom == speeds(idx_sloom);
    uncensored  = mask_s & ~bl.is_censored & ~isnan(avg_sm_corrected);
    censored    = mask_s &  bl.is_censored & ~isnan(avg_sm_corrected);

    figure; hold on
    scatter(avg_sm_corrected(uncensored), duration_s(uncensored), 8, [0.7 0.7 0.7], 'filled')
    scatter(avg_sm_corrected(censored),   duration_s(censored),   8, [0.9 0.5 0.5], 'filled', ...
        'MarkerFaceAlpha', 0.5)   % censored bouts in red — duration is a lower bound

    % Per-condition means (uncensored only for the group markers)
    for c = 1:length(conditions)
        idx_c = uncensored & bl.moving_flies == conditions(c);
        if any(idx_c)
            plot(mean(avg_sm_corrected(idx_c)), mean(duration_s(idx_c)), 'o', ...
                'MarkerFaceColor', cmap(c,:), 'MarkerEdgeColor', 'k', 'MarkerSize', 10)
        end
    end

    % Correlation on uncensored bouts only (censored durations are lower bounds)
    [r, p] = corr(avg_sm_corrected(uncensored), duration_s(uncensored), 'Type', 'Spearman');
    xlabel('Mean social motion during freeze (corrected)')
    ylabel('Freeze duration (s)')
    legend(['Uncensored', 'Censored (lower bound)', ...
        arrayfun(@(c) sprintf('%d moving', c), conditions, 'UniformOutput', false)], ...
        'Location', 'northeast')
    title(sprintf('%d cm/s  —  r = %.2f, p = %.3f  (n=%d uncensored, %d censored)', ...
        speeds(idx_sloom), r, p, sum(uncensored), sum(censored)))

    fprintf('\n%d cm/s — censored: %d / %d bouts (%.1f%%)\n', speeds(idx_sloom), ...
        sum(mask_s & bl.is_censored), sum(mask_s), ...
        100 * mean(bl.is_censored(mask_s)))
end

% ── Summary stats ─────────────────────────────────────────────────────────
fprintf('\nMedian freeze duration (s) by condition:\n')
%fprintf('  %8s', arrayfun(@(s) sprintf('%dcm/s', s), speeds, 'UniformOutput', false){:})
fprintf('\n')
for c = 1:length(conditions)
    fprintf('  %d moving:', conditions(c))
    for idx_sloom = 1:length(speeds)
        idx = bl.moving_flies == conditions(c) & bl.sloom == speeds(idx_sloom);
        fprintf('  %5.2f s ', median(duration_s(idx), 'omitnan'))
    end
    fprintf('\n')
end
