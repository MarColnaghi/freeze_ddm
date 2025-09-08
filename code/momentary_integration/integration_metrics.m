threshold_imm = 3;
threshold_mob = 3;
threshold_pc = 4;
id_code = sprintf('imm%d_mob%d_pc%d', threshold_imm, threshold_mob, threshold_pc);

exporting = false;
col = cmapper();

paths = path_generator('folder', 'integration/early_vs_late', 'bouts_id', id_code);

if isfile(fullfile(paths.dataset, 'bouts.mat'))
    s = load(fullfile(paths.dataset, 'bouts.mat'));
    bouts = s.bouts;
    s = load(fullfile(paths.dataset, 'soc_mot.mat'));
    soc_mot = s.soc_mot;

else
    disp('Could not find the bouts file. Suggest running the function generate_dataset.m another time')
end

% Set the thresholds for the metrics of interest

close all
min_dur = 0.35;
min_sm = 0.1;
max_sm = 0.7;
opt.thr_slope = 0.0;
opt.thr_com = 0.42; % Set 1 if you don't want it
opt.mov_mean = 1; % Set 1 if you don't want it
opt.thr_trend = 0.55; % Set 1 if you don't want it
opt.thr_es = .8;

for idx_slooms = [unique(bouts.sloom)]'

    if idx_slooms == 25
        sloom = 'slow';
    elseif idx_slooms == 50
        sloom = 'fast';
    end

    % Select Freezes
    bouts_proc = data_parser_new(bouts, 'type', 'immobility', 'period','loom', 'window','le');

    % Only above than minimum duration
    bouts_proc = bouts_proc(bouts_proc.durations_s >= min_dur, :);
    bouts_proc = bouts_proc(bouts_proc.sloom == idx_slooms, :);

    bouts_proc_filtered = bouts_proc(bouts_proc.avg_sm_freeze_norm >= min_sm & bouts_proc.avg_sm_freeze_norm <= max_sm, :);

    % Select freezes
    [a, ~] = ismember(soc_mot.id, bouts_proc_filtered.id);
    soc_mot_proc = soc_mot(a, :);

    % Run selection function
    [early_sm, late_sm, com, slopes, Signal_Slope, Center_of_Mass, Trend_Score, ES_difference] =  partition_variable_freezes(soc_mot_proc.ts_sm, opt);

    % All the plots
    fh = figure('color', 'w', 'Position', [100, 100, 1500, 300]);
    tl = tiledlayout(1, 5, 'TileSpacing', 'compact', 'Padding', 'loose');

    nexttile([1, 1])
    hold on
    histogram(bouts_proc.avg_sm_freeze_norm, 0:0.025:2, 'FaceColor', 'k', 'FaceAlpha', 0.1, 'EdgeColor', 'none')
    histogram(bouts_proc_filtered.avg_sm_freeze_norm, 0:0.025:2, 'FaceColor', 'k', 'FaceAlpha', 0.2, 'EdgeColor', 'none')
    histogram(bouts_proc_filtered.avg_sm_freeze_norm(late_sm), min_sm:0.025:max_sm, 'FaceColor', col.late_sm, ...
        'EdgeColor', 'none')
    histogram(bouts_proc_filtered.avg_sm_freeze_norm(early_sm), min_sm:0.025:max_sm, 'FaceColor', col.early_sm, ...
        'EdgeColor', 'none')
    xline([min_sm max_sm], '--', {num2str(min_sm), num2str(max_sm)}, 'FontSize', 16, 'LabelVerticalAlignment',  'top')
    text(1, 200, ['early - ' sprintf('n: %d', length(early_sm))], 'Color', col.early_sm, 'FontSize', 16, 'HorizontalAlignment', 'left')
    text(1, 180, ['late - ' sprintf('n: %d', length(late_sm))], 'Color', col.late_sm, 'FontSize', 16, 'HorizontalAlignment', 'left')
    ax = gca;
    apply_generic(ax, 20)
    xlabel('Avg. Social Motion')
    ylabel('Count')

    nexttile([1, 1])
    scatters(Center_of_Mass, Signal_Slope, early_sm, late_sm, opt.thr_com, opt.thr_slope, col)

    nexttile([1, 1])
    scatters(Center_of_Mass, Trend_Score, early_sm, late_sm, opt.thr_com, opt.thr_trend, col)

    nexttile([1, 1])
    scatters(Trend_Score, Signal_Slope, early_sm, late_sm, opt.thr_trend, opt.thr_slope, col)

    nexttile([1, 1])
    scatters(Center_of_Mass, ES_difference, early_sm, late_sm, opt.thr_com, opt.thr_es, col)

    exporter(fh, paths, sprintf('1_metrics_%s.pdf', sloom), 'flag', exporting)
    exporter(fh, paths, sprintf('1_metrics_%s.png', sloom), 'flag', exporting)

    fh = figure('color', 'w', 'Position', [100, 100, 800, 600]);
    tl = tiledlayout(3, 3, 'TileSpacing', 'loose', 'Padding', 'loose');

    nexttile(1, [1, 1])
    histograms(bouts_proc_filtered, 'avg_sm_freeze_norm',  early_sm, late_sm, col, 0.05)

    nexttile(4, [1, 1])
    histograms(bouts_proc_filtered, 'avg_fs_1s_norm',  early_sm, late_sm, col, 0.1)

    nexttile(7, [1, 1])
    histograms(bouts_proc_filtered, 'nloom',  early_sm, late_sm, col, 5)

    nexttile(2, [1 1])
    timeseries(soc_mot_proc, early_sm, col)

    nexttile(3, [1 1])
    timeseries(soc_mot_proc, late_sm, col)

    nexttile(5, [2 2])
    histograms(bouts_proc_filtered, 'durations_s',  early_sm, late_sm, col, 1)

    exporter(fh, paths, sprintf('2_fds_%s.pdf', sloom), 'flag', exporting)
    exporter(fh, paths, sprintf('2_fds_%s.png', sloom), 'flag', exporting)

    % If you want the imagesc as well :)

    soc_mot_proc_selected = soc_mot_proc.ts_sm(late_sm);
    vec_lengths = cellfun(@length, soc_mot_proc_selected);

    % Sort by duration
    [sorted_lengths, sort_idx] = sort(vec_lengths, 'descend');

    % Reorder the cell array accordingly
    sorted_cells = soc_mot_proc_selected(sort_idx);

    aligned_mat = padder_for_imagesc(sorted_cells, sorted_lengths, 'offset');
    fh = plot_sta_sm(aligned_mat, 'direction', 'offset', 'clim', [0 15], 'center', 'mean');

    exporter(fh, paths, sprintf('3_late_ts_%s.pdf', sloom))
    exporter(fh, paths, sprintf('3_late_ts_%s.png', sloom))

    % for early trials
    soc_mot_proc_selected = soc_mot_proc.ts_sm(early_sm);
    vec_lengths = cellfun(@length, soc_mot_proc_selected);

    % Sort by duration
    [sorted_lengths, sort_idx] = sort(vec_lengths, 'descend');

    % Reorder the cell array accordingly
    sorted_cells = soc_mot_proc_selected(sort_idx);

    aligned_mat = padder_for_imagesc(sorted_cells, sorted_lengths, 'onset');
    fh = plot_sta_sm(aligned_mat, 'direction', 'onset', 'clim', [0 15], 'center', 'mean');
    ylim([0 max([length(late_sm), length(early_sm)])])

    exporter(fh, paths, sprintf('3_early_ts_%s.pdf', sloom), 'flag', exporting)
    exporter(fh, paths, sprintf('3_early_ts_%s.png', sloom), 'flag', exporting)

end

% Succesful Combinations:

% % # 1 - 16.01.2025
% min_dur = 0.35;
% min_sm = 0.1;
% max_sm = 0.8;
% opt.thr_slope = 0.0;
% opt.thr_com = 0.4; % Set 1 if you don't want it
% opt.mov_mean = 6; % Set 1 if you don't want it
% opt.thr_trend = 0.52; % Set 1 if you don't want it
% opt.thr_es = 3;
%
% % # 2 - 16.01.2025
% min_dur = 0.35;
% min_sm = 0.1;
% max_sm = 0.6;
% opt.thr_slope = 0.0;
% opt.thr_com = 0.4; % Set 1 if you don't want it
% opt.mov_mean = 1; % Set 1 if you don't want it
% opt.thr_trend = 0.52; % Set 1 if you don't want it
% opt.thr_es = -20;
%
% % # 3 - 16.01.2025 - I really like this - Speed Division as well :)
% min_dur = 0.35;
% min_sm = 0.1;
% max_sm = 0.6;
% opt.thr_slope = 0.0;
% opt.thr_com = 0.45; % Set 1 if you don't want it
% opt.mov_mean = 1; % Set 1 if you don't want it
% opt.thr_trend = 0.52; % Set 1 if you don't want it
% opt.thr_es = 1;
%
% % # 4 - 25.01.2025
% min_dur = 0.35;
% min_sm = 0.1;
% max_sm = 0.7;
% opt.thr_slope = 0.0;
% opt.thr_com = 0.42; % Set 1 if you don't want it
% opt.mov_mean = 1; % Set 1 if you don't want it
% opt.thr_trend = 0.55; % Set 1 if you don't want it
% opt.thr_es = .8;


function scatters(xdata, ydata, early_idx, late_idx, x_thr, y_thr, col)

scatter(xdata, ydata, 8, 'filled', 'o', ...
    'MarkerFaceColor', 'k', 'MarkerFaceAlpha', 0.15, 'MarkerEdgeColor', 'none')
hold on
scatter(xdata(late_idx), ydata(late_idx), 8, 'filled', 'o', ...
    'MarkerFaceColor', col.late_sm , 'MarkerFaceAlpha', 0.5, 'MarkerEdgeColor', 'none')
scatter(xdata(early_idx), ydata(early_idx), 8, 'filled', 'o', ...
    'MarkerFaceColor', col.early_sm , 'MarkerFaceAlpha', 0.5, 'MarkerEdgeColor', 'none')

xline([x_thr, 1 - x_thr], '--', {num2str(x_thr), num2str(1 - x_thr)}, ...
    'FontSize', 16, 'LabelVerticalAlignment', 'bottom')

if strcmp(inputname(2), 'Trend_Score')
    yline([1 - y_thr, y_thr], '--', {num2str(y_thr), ' '}, 'FontSize', 16)
    ylim([0 1])
else
    yline([-y_thr, y_thr], '--', {num2str(y_thr), ' '}, 'FontSize', 16)
    ylim([-round(max(ydata) * 1.3) round(max(ydata) * 1.3)])
end

xlim([0 1])
ax = gca;
apply_generic(ax, 20)

xlabel(strrep(inputname(1), '_', ' '))
ylabel(strrep(inputname(2), '_', ' '))

end


function histograms(bouts_proc_filtered, variable, early_idx, late_idx, col, bin_size)

hold on
data = bouts_proc_filtered.(variable)([early_idx late_idx]);
mins = min(data); maxs = max(data) + bin_size;

data = bouts_proc_filtered.(variable)(early_idx);
histogram(data, mins:bin_size:maxs, 'FaceColor', col.early_sm, ...
    'EdgeColor', 'none', 'Normalization', 'probability', 'FaceAlpha', 0.7)

data = bouts_proc_filtered.(variable)(late_idx);
histogram(data, mins:bin_size:maxs, 'FaceColor', col.late_sm, ...
    'EdgeColor', 'none', 'Normalization', 'probability', 'FaceAlpha', 0.7)

xlim([mins - 0.1 maxs + 0.1])
ax = gca;
xticks([round(mins, 1) round(maxs,1)])
ylabel('Fraction')

if strcmp(variable, 'avg_sm_freeze_norm')
    xlabel('Social Motion')
    ylim([-ax.YAxis.Limits(2)/20 ax.YAxis.Limits(2)])

elseif strcmp(variable, 'avg_fs_1s_norm')
    xlabel('Focal Speed')
    ylim([-ax.YAxis.Limits(2)/20 ax.YAxis.Limits(2)])

elseif strcmp(variable, 'nloom')
    xlabel('Loom Number')
    xticks(3:5:18)
    xticklabels({'1-5', '6-10', '11-15', '16-20'})
    xtickangle(30)
    xlim([0 22])
    ylim([-ax.YAxis.Limits(2)/20 ax.YAxis.Limits(2)])

elseif strcmp(variable, 'durations_s')
    xlabel('Freeze Duration (s)')
    xlim([0 31])
    xticks([0 10 20 30])
    ylim([-ax.YAxis.Limits(2)/40 ax.YAxis.Limits(2)])
    ylim([-0.01 0.55])

end
apply_generic(ax, 20)

end

function timeseries(soc_mot_proc, idx, col)
hold on
ax = gca;
ax.YAxis.Visible = 'off';
xticks([0 600 1200])

if strcmp(inputname(2), 'early_sm')
    cellfun(@(x) plot(1:length(x), x, 'Color', [hex2rgb(col.early_sm) 0.08]), soc_mot_proc.ts_sm(idx));
    xticks([0 600 1200])
    xticklabels({'Freeze Entry', '10', '20'})
    xtickangle(0)
    ylim([-3 63])
    xlim([0 600])
    xlabel('Time from Entry (s)')

elseif strcmp(inputname(2), 'late_sm')
    cellfun(@(x) plot(-length(x):-1, x, 'Color', [hex2rgb(col.late_sm) 0.08]), soc_mot_proc.ts_sm(idx));
    xticks([-1200 -600 -1])
    xticklabels({'-20', '-10', 'Freeze Exit'})
    xtickangle(0)
    ylim([-3 63])
    xlim([-600 0])
    xlabel('Time from Exit (s)')

end

apply_generic(ax, 20)

end
