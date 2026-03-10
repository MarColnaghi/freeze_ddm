% --- 1. Setup and Data Loading ---
model_list = {'dddim2'};
run_list   = {'run01_010326_1800'};

% Adjust these paths to match your local PhD directory structure if needed
version = {'_v8'};

paths_analysis = path_generator('folder', fullfile('momentary_integration','alternative_test', model_list{1}, strcat(run_list{1}, version{1})));
paths_analysis.fig = fullfile(paths_analysis.fig, 'detailed_analysis');
col = cmapper('', 2);

% Load Likelihoods and Data
ll = struct2table(importdata(fullfile(paths_analysis.results, 'll_table.mat')));
data = alternative_test_collect_data(model_list, run_list);

% Load Motion Cache (from your analyse_sm logic)
% Assuming data(1).sim_params contains the path to motion_cache
motion_cache = importdata(fullfile(paths_analysis.cache_path, 'motion_cache.mat'));

% --- 2. Configuration ---
align_to = 'onset'; % 'onset' or 'offset'
trace_len = numel(data.t);    % Frames to show

freezes_ref = data(1).table;

% Remove censored freezes
is_censored = freezes_ref.durations_s > data(1).results.points.censoring;
freezes_ref = freezes_ref(~ is_censored,:);
ll = ll(~ is_censored,:);
data(1).table = data(1).table(~ is_censored,:);
data(1).extra_tv = data(1).extra_tv(~ is_censored,:);

n_bouts     = height(freezes_ref);
diff = ll.tv - ll.st;
[sorted_deltall, sort_idx] = sort(diff);

value_xlines = 0.25;

idx_crossing(2) = find(sorted_deltall > 0, 1);
idx_crossing(1) = find(sorted_deltall > -value_xlines, 1);
idx_crossing(3) = find(sorted_deltall > value_xlines, 1);

hold on
ylim([-2 2])
xlim([-50 height(ll) + 50])
ax = gca;

opt.thr_slope = 0.0;
opt.thr_com = 0.42; % Set 1 if you don't want it
opt.mov_mean = 1; % Set 1 if you don't want it
opt.thr_trend = 0.55; % Set 1 if you don't want it
opt.thr_es = .8;

cell_sm = extract_sm_from_bouts(freezes_ref, 'formula', 'cell', 'type', 'onlyfreeze');
[early_sm, late_sm, com, slopes, Signal_Slope, Center_of_Mass, Trend_Score, ES_difference, Variance] =  partition_variable_freezes(cell_sm, opt);

metrics = {'Signal_Slope', 'Center_of_Mass', 'ES_difference', 'Variance'};
limits = {[-0.1 0.1], [0 1], [-2 2], [-0.05 1.05], [0.45 0.55]};
xticks = {[-0.1 0 0.1], [0 0.5 1], [-2 0 2], [0 0.5 1], [0.45 0 0.55]};

% --- 3. Figure and Plotting Loop ---
% Ensure sort_ll is defined before the loop
sort_ll = diff(sort_idx);

marker_size = 10;

for idx_metrics = 1:length(metrics)

    % 1. Extract the metric and sort it based on Delta-LL
    % Since partition_variable_freezes was run on the unsorted cell_sm,
    % we must apply the sort_idx here.
    current_metric_name = metrics{idx_metrics};
    raw_metric_data = eval(current_metric_name); % Get the variable by string name
    sorted_metric = raw_metric_data(sort_idx);
    sorted_durations = freezes_ref.durations_s(sort_idx);

    % 2. Initialize Figure
    fh = figure('color', 'w', 'Position', [100, 100, 1100, 350]);
    tl = tiledlayout(1, 3, 'TileSpacing', 'loose', 'Padding', 'compact');

    % Panel A: Metric vs. Rank (Sorted by Delta-LL)
    ax(1) = nexttile;
    scatter(1:n_bouts, sorted_metric, marker_size, sorted_durations, 'filled', 'MarkerFaceAlpha', 0.4);
    clim([0 5]);
    ylabel(strrep(current_metric_name, '_', ' '));
    xlabel('Rank (st > tv)');
    apply_generic(ax(1), 'font_size', 18, 'ylims', limits{idx_metrics}, 'xlims', [-50 n_bouts + 50]);
    xline(idx_crossing(1), 'k--')
    xline(idx_crossing(2), 'k-.')
    xline(idx_crossing(3), 'k-.')

    % Panel B: Metric vs. Delta-LL (Black)
    ax(2) = nexttile;
    scatter(sort_ll, sorted_metric, marker_size, 'k', 'filled', 'MarkerFaceAlpha', 0.4);
    xlabel('$\log \!\big(\mathcal{L}_{\mathrm{tv}} / \mathcal{L}_{\mathrm{st}}\big)$', 'Interpreter','latex');
    yline(sum(limits{idx_metrics})/2, '--', 'Alpha', 0.5);
    xline(-value_xlines, 'k--')
    xline(0, 'k-.')
    xline(value_xlines, 'k-.')
    apply_generic(ax(2), 'font_size', 18, 'ylims', limits{idx_metrics}, 'xlims', [-2, 2], 'no_yticks', true);

    % Panel C: Metric vs. Delta-LL (Colored by Duration)
    ax(3) = nexttile;
    scatter(sort_ll, sorted_metric, marker_size, sorted_durations, 'filled', 'MarkerFaceAlpha', 0.4);
    clim([0 5]);
    colorcet('I1'); % Requires colorcet toolbox
    xlabel('$\log \!\big(\mathcal{L}_{\mathrm{tv}} / \mathcal{L}_{\mathrm{st}}\big)$', 'Interpreter','latex');
    yline(sum(limits{idx_metrics})/2, '--', 'Alpha', 0.5);
    xline(-value_xlines, 'k--')
    xline(0, 'k-.')
    xline(value_xlines, 'k-.')
    apply_generic(ax(3), 'font_size', 18, 'ylims', limits{idx_metrics}, 'xlims', [-2, 2], 'no_ylabels', true);

    % Add Colorbar to the final panel
    cb = colorbar;
    cb.Label.String = 'Freeze Duration (s)';
    cb.LineWidth = 2;
    cb.FontSize = 18;

    linkaxes(ax(:), 'y')
    exporter(fh, paths_analysis, sprintf('%s.pdf', metrics{idx_metrics}));

end