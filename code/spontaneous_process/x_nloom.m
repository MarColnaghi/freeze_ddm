id_code = 'imm2_mob2_pc4';
col = cmapper();
paths_out = path_generator('folder', '/spontaneous_process', 'bouts_id', id_code, 'imfirst', false);
bouts = importdata(fullfile(paths_out.dataset, 'bouts.mat'));
bouts_spontaneous = data_parser_new(bouts, 'period', 'bsl', 'window', 'le', 'type', 'immobility', 'nloom', 10:20);
bouts_le = data_parser_new(bouts, 'period', 'loom', 'window', 'le', 'type', 'immobility', 'nloom', 2:20);

col = cmapper('', 30);
fh = figure('Position', [100 100 400 400], 'Color', 'w');
t = tiledlayout(1, 1, 'TileSpacing', 'compact', 'Padding', 'loose');
hold on

for idx_loom = 1:20
    bouts_spontaneous = data_parser_new(bouts, 'period', 'bsl', 'window', 'all', 'type', 'immobility', 'nloom', idx_loom);
    [f, x] = ecdf(bouts_spontaneous.durations); 
    plot(x, f, 'Color', col.vars.nloom(10 + idx_loom,:), 'LineWidth', 1);
end

apply_generic(gca, 'xlim', [0 60], 'ylim', [-0.05 1.05])
xlabel('Immobility Frames')
ylabel('ecdf')

exporter(fh, paths_out, 'contaminant_nloom.pdf')