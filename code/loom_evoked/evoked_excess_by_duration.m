% Loom-evoked excess against freeze duration. For each range of freeze duration,
% the onsets that fall inside the window under test (loom_evoked_data) are counted
% for the real looms and for the baseline mock looms. A ratio above 1 is what the
% loom added, and the first range from which it stays above 1 says how short a
% freeze can be and still belong to the response, which is what the minimum
% duration of loom_evoked_bouts assumes. The dashed line is the control: onsets
% long after the loom, where there is no response left, so it shows the ratio the
% two periods have anyway.
%
% No length rule and no contact censoring are applied here, since the length rule
% is what is being tested, so the durations are the real ones. Flies that freeze a
% lot at baseline are excluded, as in loom_evoked_bouts.

[bouts, thresholds] = loom_evoked_data;

% Freeze duration ranges in frames, both ends included, cut at these breakpoints.
% One of them is the minimum duration of loom_evoked_bouts, so a boundary always
% sits on it
rule_min = thresholds.min_dur;
breaks = unique([0 6 8 10 rule_min 14 16 19 23 30 45 60 120 300]);
duration_ranges = [breaks(1:end-1)', breaks(2:end)' - 1; breaks(end), Inf];

% Onsets this far after the loom are the control
far_latency = [150 540];

bouts_all = data_parser_new(bouts, 'type', 'immobility', 'min_dur', 0, 'exclude_flies', true);

speeds = [25 50];
windows = {thresholds.le_window_sl, thresholds.le_window_fl};
n_ranges = size(duration_ranges, 1);

n_loom = zeros(n_ranges, 2);
n_base = zeros(n_ranges, 2);
far_loom = zeros(n_ranges, 2);
far_base = zeros(n_ranges, 2);

is_loom = bouts_all.period == 1;
is_base = bouts_all.period == 0;
in_far = bouts_all.onsets_loomwin >= far_latency(1) & bouts_all.onsets_loomwin <= far_latency(2);

for i = 1:2
    in_speed = bouts_all.sloom == speeds(i);
    in_window = bouts_all.onsets_loomwin >= windows{i}(1) & bouts_all.onsets_loomwin <= windows{i}(2);

    for r = 1:n_ranges
        in_range = bouts_all.durations >= duration_ranges(r, 1) & bouts_all.durations <= duration_ranges(r, 2);

        n_loom(r, i) = sum(in_speed & in_range & in_window & is_loom);
        n_base(r, i) = sum(in_speed & in_range & in_window & is_base);
        far_loom(r, i) = sum(in_speed & in_range & in_far & is_loom);
        far_base(r, i) = sum(in_speed & in_range & in_far & is_base);
    end
end

% Ratio with a 95% interval from the Poisson error of the two counts. Half a
% count is added to both, so that a range with no baseline onsets still gives a
% finite ratio and interval
ratio = (n_loom + 0.5) ./ (n_base + 0.5);
log_se = sqrt(1 ./ (n_loom + 0.5) + 1 ./ (n_base + 0.5));
ci_low = ratio .* exp(-1.96 * log_se);
ci_high = ratio .* exp(1.96 * log_se);
excess = n_loom - n_base;
far_ratio = far_loom ./ far_base;

range_label = cell(n_ranges, 1);
for r = 1:n_ranges
    if isinf(duration_ranges(r, 2))
        range_label{r} = sprintf('%d+', duration_ranges(r, 1));
    else
        range_label{r} = sprintf('%d-%d', duration_ranges(r, 1), duration_ranges(r, 2));
    end
end

for i = 1:2
    fprintf('\n%d cm/s, window [%d %d]\n', speeds(i), windows{i}(1), windows{i}(2));
    disp(table(range_label, n_loom(:, i), n_base(:, i), ratio(:, i), ci_low(:, i), ci_high(:, i), excess(:, i), far_ratio(:, i), ...
        'VariableNames', {'duration', 'loom', 'baseline', 'ratio', 'ci_low', 'ci_high', 'excess', 'far_ratio'}))

    % First range from which the lower end of the interval stays above 1
    first = find(~(ci_low(:, i) > 1), 1, 'last') + 1;
    if isempty(first), first = 1; end
    if first > n_ranges
        fprintf('no range keeps its interval above 1 to the end\n');
    else
        fprintf('the interval stays above 1 from the %s frame range on\n', range_label{first});
    end
end

col = cmapper([], 5);
fh = figure('color', 'w', 'Position', [100, 100, 900, 750]);
tl = tiledlayout(fh, 2, 1, 'TileSpacing', 'compact');
x_rule = find(duration_ranges(:, 1) == rule_min) - 0.5;

ax(1) = nexttile;
hold on
h = gobjects(1, 4);
for i = 1:2
    c = col.vars.sloom(3 * i, :);
    h(2*i - 1) = errorbar(1:n_ranges, ratio(:, i), ratio(:, i) - ci_low(:, i), ci_high(:, i) - ratio(:, i), '-o', 'Color', c, 'MarkerFaceColor', c, 'LineWidth', 1.5);
    h(2*i) = plot(1:n_ranges, far_ratio(:, i), '--', 'Color', c, 'LineWidth', 1.2);
end
yline(1, 'k:');
xline(x_rule, 'k-', 'minimum duration', 'LabelVerticalAlignment', 'bottom', 'FontSize', 14);
set(ax(1), 'YScale', 'log');
ylabel('loom / baseline')
legend(h, {'25 cm/s, in window', '25 cm/s, far (control)', '50 cm/s, in window', '50 cm/s, far (control)'}, 'Location', 'northwest', 'Box', 'off', 'FontSize', 14);

ax(2) = nexttile;
bars = bar(1:n_ranges, excess);
for i = 1:2
    bars(i).FaceColor = col.vars.sloom(3 * i, :);
end
hold on
yline(0, 'k-');
xline(x_rule, 'k-');
ylabel('onsets above baseline')
xlabel('Freeze duration (frames)')

for k = 1:2
    xlim(ax(k), [0.5, n_ranges + 0.5]);
    xticks(ax(k), 1:n_ranges);
    xticklabels(ax(k), range_label);
    apply_generic(ax(k), 'font_size', 16);
end
xtickangle(ax(1), 45);
xtickangle(ax(2), 45);
