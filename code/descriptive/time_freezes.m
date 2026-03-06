id_code = 'imm3_mob3_pc4';
paths_out = path_generator('folder', 'descriptive/fd_durs', 'bouts_id', id_code, 'imfirst', false);

col = cmapper();
thresholds = define_thresholds;
thresholds.le_window_fl = [5 40];
thresholds.le_window_sl = [15 50];

bouts = importdata(fullfile(paths_out.dataset, 'bouts.mat'));
bouts = bouts_formatting(bouts, thresholds);

paths_out = path_generator('folder', '/spontaneous_process', 'bouts_id', id_code, 'imfirst', false);
bouts_spontaneous = data_parser_new(bouts, 'period', 'bsl', 'window', 'le', 'type', 'immobility', 'nloom', 10:20);
bouts_le = data_parser_new(bouts, 'period', 'loom', 'window', 'le', 'type', 'immobility', 'nloom', 2:20, 'min_dur', 30);
bouts_all = data_parser_new(bouts, 'min_dur', 30);

[G, fly_ids, periods] = findgroups(bouts_all.fly, bouts_all.period);

total_freeze_time = splitapply(@sum, bouts_all.durations_s, G);

T_freeze = table(fly_ids, periods, total_freeze_time, ...
    'VariableNames', {'fly','period','total_freeze_time_s'});

T_freeze = sortrows(T_freeze, {'period','total_freeze_time_s'}, ...
                               {'descend','descend'});

T_freeze.period = categorical(T_freeze.period, ...
    [0 1], {'bsl','loom'});

T_wide = unstack(T_freeze, 'total_freeze_time_s', 'period');

threshold_duration = 100;
long_freezes_during_bsl = T_wide.bsl > threshold_duration;
sortrows(T_wide, 'bsl')

sprintf('The flies with more than %ds of freeze are: %s', threshold_duration, num2str(T_wide.fly(long_freezes_during_bsl)'))
fh = figure('Position', [100 100 500 600], 'Color', 'w'); hold on
scatter(T_wide.bsl, T_wide.loom, 24, long_freezes_during_bsl, 'filled', 'MarkerFaceAlpha', .55)
colormap(cbrewer2('RdBu', []))
clim([-.15 1.15])
xline(threshold_duration, 'r--')
xlabel('Total freeze time — Baseline (s)')
ylabel('Total freeze time — Loom (s)')
axis square
grid on
apply_generic(gca)