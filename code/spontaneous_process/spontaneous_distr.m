id_code = 'imm2_mob2_pc4';
paths_out = path_generator('folder', '/spontaneous', 'bouts_id', id_code, 'imfirst', false);
bouts = importdata(fullfile(paths_out.dataset, 'bouts.mat'));
bouts_spontaneous = data_parser_new(bouts, 'period', 'bsl', 'window', 'le', 'type', 'immobility');
bouts_le = data_parser_new(bouts, 'period', 'loom', 'window', 'le', 'type', 'immobility');

fh = figure('color', 'w', 'Position', [100, 100, 900, 350]);
hold on
histogram(bouts_spontaneous.durations, 'BinMethod', 'integers', 'Normalization', 'probability','EdgeColor','none')
histogram(bouts_le.durations, 'BinMethod', 'integers', 'Normalization', 'probability', 'EdgeColor','none')
xlim([-2 302])
ylim([-0.005 0.205])
apply_generic(gca)

fh = figure('color', 'w', 'Position', [100, 100, 900, 350]);
hold on
histogram(bouts_spontaneous.avg_sm_freeze_norm, 0:0.01:2, 'Normalization', 'probability', 'EdgeColor','none')
histogram(bouts_le.avg_sm_freeze_norm, 0:0.01:2, 'Normalization', 'probability', 'EdgeColor','none')
apply_generic(gca)

fh = figure('color', 'w', 'Position', [100, 100, 900, 350]);
hold on
histogram(bouts_spontaneous.avg_fs_1s_norm, 0:0.01:2, 'Normalization', 'probability', 'EdgeColor','none')
histogram(bouts_le.avg_fs_1s_norm, 0:0.01:2, 'Normalization', 'probability', 'EdgeColor','none')
apply_generic(gca)

col = cbrewer2('Purples', 20);
fh = figure('color', 'w', 'Position', [100, 100, 450, 450]);
hold on
for idx_nloom = 1:20
    bouts_spontaneous = data_parser_new(bouts, 'period', 'bsl', 'window', 'all', 'type', 'immobility', 'nloom', idx_nloom);
    [cdf, fd] = ecdf(bouts_spontaneous.durations);
    plot(fd, cdf, 'Color', col(idx_nloom, :), 'LineWidth', 1.2);
   
end
xlim([0 60])
ylim([0 1])
apply_generic(gca)

%%