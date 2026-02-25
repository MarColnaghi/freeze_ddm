fs = 60;
id_code = 'imm2_mob2_pc4';
col = cmapper();
paths_out = path_generator('folder', 'spontaneous_process', 'bouts_id', id_code, 'imfirst', false);

bouts = importdata(fullfile(paths_out.dataset, 'bouts.mat'));
bouts_spontaneous = data_parser_new(bouts, 'period', 'bsl', 'window', 'all', 'type', 'immobility', 'nloom', 2:20);
bouts_le = data_parser_new(bouts, 'period', 'loom', 'window', 'le', 'type', 'immobility', 'nloom', 2:20, 'min_dur', 30);

sm_spontaneous = extract_sm_from_bouts(bouts_spontaneous, 'type', 'onlyfreeze', 'align', 'offset');
sm_le = extract_sm_from_bouts(bouts_le, 'type', 'onlyfreeze', 'align', 'offset');
t_sp = (-(size(sm_spontaneous, 2)-1) : 0 )/ fs;
t_le = (-(size(sm_le, 2)-1) : 0 )/ fs;

fh = figure('color','w','Position',[100, 100, 520, 270]);
tiledlayout(1,1, 'TileSpacing', 'compact', 'Padding', 'compact')
hold on

m_sp = mean(sm_spontaneous, 'omitnan');
n_sp = sum(~isnan(sm_spontaneous), 1);
sem_sp = std(sm_spontaneous, 'omitnan') ./ sqrt(n_sp);

m_le = mean(sm_le, 'omitnan');
n_le = sum(~isnan(sm_le), 1);
sem_le = std(sm_le, 'omitnan') ./ sqrt(n_le);

hold on
fill([t_sp fliplr(t_sp)], ...
     [m_sp-sem_sp fliplr(m_sp+sem_sp)], ...
     hex2rgb(col.period.bsl), 'FaceAlpha', 0.25, 'EdgeColor', 'none')

fill([t_le fliplr(t_le)], ...
     [m_le-sem_le fliplr(m_le+sem_le)], ...
     hex2rgb(col.period.loom), 'FaceAlpha', 0.25, 'EdgeColor', 'none')

plot(t_sp, m_sp, 'LineWidth', 2, 'Color', col.period.bsl)
plot(t_le, m_le, 'LineWidth', 2, 'Color', col.period.loom)

plot([-median(bouts_spontaneous.durations_s) -median(bouts_spontaneous.durations_s)] ,[0 1], 'Color', col.period.bsl, 'LineStyle', '-.')
plot([-median(bouts_le.durations_s) -median(bouts_le.durations_s)] ,[0 1], 'Color', col.period.loom, 'LineStyle', '-.')

apply_generic(gca, 'xlim', [-10.1 0.1], 'ylim', [0 1], 'ypos', 'right', 'xticks', -10:2:0)
ylabel('Social Motion')
xlabel('Time')

exporter(fh, paths_out, 'sm_contaminant.pdf')