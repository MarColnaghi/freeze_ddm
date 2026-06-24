clearvars; close all; clc;

% Load colors
num_quantiles = 5;
extra.quantiles = num_quantiles;
col       = cmapper([], extra.quantiles);
col_nloom = cmapper([], 30);

% Load paths
paths = path_generator('folder', fullfile('descriptive', 'rasters'));

% Load cached timeseries
sm_cache   = importdata(fullfile(paths.cache_path, 'motion_cache.mat'));
pc_cache   = importdata(fullfile(paths.cache_path, 'pixel_cache.mat'));
fr_cache   = importdata(fullfile(paths.cache_path, 'freeze_cache.mat'));
fs_cache   = importdata(fullfile(paths.cache_path, 'speed_cache.mat'));
md_cache   = importdata(fullfile(paths.cache_path, 'mindist_cache.mat'));
loom_cache = importdata(fullfile(paths.cache_path, 'loom_cache.mat'));

% Load bouts
threshold_imm = 2; threshold_mob = 2; threshold_pc = 4;
id_code    = sprintf('imm%d_mob%d_pc%d', threshold_imm, threshold_mob, threshold_pc);
thresholds = define_thresholds;
bouts      = importdata(fullfile(paths.dataset, 'bouts.mat'));
bouts      = bouts_formatting(bouts, thresholds);

% Select flies at the chosen loom speed
ls             = 25;
selected_flies = unique(bouts.fly(bouts.sloom == ls, :));
n_moving_flies = accumarray(bouts.fly, bouts.moving_flies, [], @(x) x(1));
n_moving_flies = n_moving_flies(selected_flies);

% Extract matrices for selected flies
sm_mat   = cache2mat(sm_cache,   'selected_flies', selected_flies');
fs_mat   = cache2mat(fs_cache,   'selected_flies', selected_flies');
fr_mat   = cache2mat(fr_cache,   'selected_flies', selected_flies');
pc_mat   = cache2mat(pc_cache,   'selected_flies', selected_flies');
md_mat   = cache2mat(md_cache,   'selected_flies', selected_flies');
loom_mat = cache2mat(loom_cache, 'selected_flies', selected_flies');

% Build summary table and sort
% n_moving_flies and md_mat rows share the same order (selected_flies).
% ts is sorted for display; use n_moving_flies (not ts.moving_flies) to
% index md_mat in analyses below.
ts              = table();
ts.freeze_time  = sum(fr_mat(:, 18000:end), 2);
ts.moving_flies = n_moving_flies;
ts.freeze       = ~fr_mat;
ts.sm           = sm_mat;
ts.fly          = selected_flies;
ts = sortrows(ts, {'moving_flies', 'freeze_time'}, 'ascend');

% Loom onset times
loom_ts         = diff(loom_mat, [], 2) == 1;
n_flies         = size(loom_mat, 1);
n_looms_per_fly = sum(loom_ts, 2);
n_looms         = n_looms_per_fly(1);
assert(all(n_looms_per_fly == n_looms), 'Flies do not all have the same number of looms.');

loom_times = nan(n_flies, n_looms);
for f = 1:n_flies
    loom_times(f, :) = find(loom_ts(f, :));
end
loom_times = [loom_times, 36000 * ones(n_flies, 1)];   % sentinel: end-of-recording frame
nIntervals = n_looms;

% Min-distance slices in a window after each loom onset
window    = 800;
md_slices = arrayfun(@(fly, k) ...
    md_mat(fly, loom_times(fly,k) : min(loom_times(fly,k)+window, size(md_mat,2))), ...
    repmat((1:n_flies)', 1, nIntervals), ...
    repmat(1:nIntervals,  n_flies,    1), ...
    'UniformOutput', false);

% Contact detection
threshold_contact = 70;
contact = cellfun(@(x) any(x < threshold_contact), md_slices);

figure
plot(sum(contact, 2))
xlabel('Fly index')
ylabel('Number of contacts')
title('Contact events per fly')

% First-contact latency
first_contact_frame          = NaN(size(md_slices));
first_contact_frame(contact) = cellfun( ...
    @(x) find(x < threshold_contact, 1, 'first'), ...
    md_slices(contact));

latency_seconds = first_contact_frame / 60;
mean_latency    = mean(latency_seconds, 1, 'omitnan');

% Control: same procedure on random baseline windows (frames 1–18000)
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

edges = (0:1/60:10.5) - 1/120;
figure; hold on
histogram(latency_seconds(:),   edges, 'Normalization', 'probability', ...
    'FaceColor', [0.2 0.5 0.8], 'FaceAlpha', 0.6, 'EdgeColor', 'none')
histogram(latency_baseline_s(:), edges, 'Normalization', 'probability', ...
    'FaceColor', [0.7 0.7 0.7], 'FaceAlpha', 0.6, 'EdgeColor', 'none')
xlabel('Latency to contact (s)')
ylabel('Probability')
legend('Loom-triggered', 'Baseline (random windows)')
fprintf('Baseline contact rate:     %.1f%%\n', 100*mean(~isnan(latency_baseline_s(:))))
fprintf('Loom contact rate:         %.1f%%\n', 100*mean(~isnan(latency_seconds(:))))

% Distance at and just before loom onset
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

% Loom-triggered average: ALL (fly, loom) pairs
% Pre-onset window is critical — if distance is already decreasing before t=0,
% the dip is NOT caused by the loom (flies were already approaching).
pre  = 180;   % 2 s before loom onset
post = 630;   % 5 s after loom onset
t_axis  = -pre:post;
n_pairs = n_flies * nIntervals;
traj_all = NaN(n_pairs, pre+post+1);
pair = 0;
for f = 1:n_flies
    for k = 1:nIntervals
        pair      = pair + 1;
        t0        = loom_times(f, k);
        t_start   = max(t0 - pre,  1);
        t_end     = min(t0 + post, size(md_mat, 2));
        rel_start = t_start - t0 + pre + 1;
        rel_end   = t_end   - t0 + pre + 1;
        traj_all(pair, rel_start:rel_end) = md_mat(f, t_start:t_end);
    end
end

mu  = mean(traj_all, 1, 'omitnan');
sem = std(traj_all,  0, 1, 'omitnan') ./ sqrt(sum(~isnan(traj_all), 1));

figure; hold on
fill([t_axis, fliplr(t_axis)], [mu+sem, fliplr(mu-sem)], [0.7 0.7 0.7], 'EdgeColor', 'none')
plot(t_axis, mu, 'k', 'LineWidth', 1.5)
xline(0, '--r', 'LineWidth', 1)
xlabel('Frames relative to loom onset')
ylabel('Mean min inter-fly distance')

% Contact rate by condition (number of moving conspecifics)
conditions    = unique(n_moving_flies)';          % e.g. [0 1 2 3 4]
contact_rate  = mean(contact, 2);                 % fraction of looms with a contact, per fly

group_means = NaN(1, length(conditions));
group_sems  = NaN(1, length(conditions));
group_n     = NaN(1, length(conditions));
for c = 1:length(conditions)
    idx             = n_moving_flies == conditions(c);
    rates           = contact_rate(idx);
    group_means(c)  = mean(rates);
    group_sems(c)   = std(rates) / sqrt(sum(idx));
    group_n(c)      = sum(idx);
end

figure; hold on
for c = 1:length(conditions)
    idx = n_moving_flies == conditions(c);
    scatter(conditions(c) * ones(group_n(c), 1) + 0.05*randn(group_n(c),1), ...
            contact_rate(idx), 15, [0.7 0.7 0.7], 'filled')
end
errorbar(conditions, group_means, group_sems, 'ko-', 'LineWidth', 1.5, 'MarkerFaceColor', 'k')
xticks(conditions)
xlabel('Number of moving conspecifics')
ylabel('Contact rate (fraction of looms with contact)')
title(sprintf('threshold = %d px, window = %.1f s', threshold_contact, window/60))

% Loom-triggered average per condition
figure; hold on
cmap = lines(length(conditions));
for c = 1:length(conditions)
    idx = find(n_moving_flies == conditions(c));
    pairs_c = [];
    for fi = 1:length(idx)
        f = idx(fi);
        for k = 1:nIntervals
            t0        = loom_times(f, k);
            t_start   = max(t0 - pre,  1);
            t_end     = min(t0 + post, size(md_mat, 2));
            row       = NaN(1, pre+post+1);
            rel_start = t_start - t0 + pre + 1;
            rel_end   = t_end   - t0 + pre + 1;
            row(rel_start:rel_end) = md_mat(f, t_start:t_end);
            pairs_c   = [pairs_c; row]; %#ok<AGROW>
        end
    end
    mu_c = mean(pairs_c, 1, 'omitnan');
    plot(t_axis, mu_c, 'Color', cmap(c,:), 'LineWidth', 1.5)
end
xline(0, '--k', 'LineWidth', 1)
yline(threshold_contact, ':k', 'LineWidth', 1)
xlabel('Frames relative to loom onset')
ylabel('Mean min inter-fly distance')
legend(arrayfun(@(c) sprintf('%d moving', c), conditions, 'UniformOutput', false))

% Median distance over time by condition
median_loom_ts = median(loom_times(:, 1:nIntervals), 1);

figure; hold on
for c = 1:length(conditions)
    plot(median(md_mat(n_moving_flies == conditions(c), :), 1), 'Color', cmap(c,:))
end
xline(median_loom_ts)
xlabel('Frame')
ylabel('Median min inter-fly distance')
legend(arrayfun(@(c) sprintf('%d moving', c), conditions, 'UniformOutput', false))
