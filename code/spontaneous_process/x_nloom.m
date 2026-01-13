% This script plots the ecdf of the spontaneous process as a function of
% loom number

id_code = 'imm2_mob2_pc4';
paths_out = path_generator('folder', '/spontaneous_process', 'bouts_id', id_code, 'imfirst', false);
bouts = importdata(fullfile(paths_out.dataset, 'bouts.mat'));

col = cmapper('', 30);
fh = figure('Position', [100 100 500 600], 'Color', 'w');
t = tiledlayout(3, 1, 'TileSpacing', 'compact', 'Padding', 'loose');
bouts_spontaneous = data_parser_new(bouts, 'period', 'bsl', 'window', 'all', 'type', 'immobility');

hold on
nexttile(1, [1,1])
bouts_spontaneous = bouts_spontaneous(~isnan(bouts_spontaneous.onsets_loomwin), :);
bh = bar(histcounts(bouts_spontaneous.nloom), 'FaceColor', 'flat', 'EdgeColor', 'none', 'CData', col.vars.nloom(11:30,:));
apply_generic(gca, 'xticks', [1:5:20 20], 'font_size', 20)
ylabel('Count')
xlabel({'Baseline Epochs', '(surrogate looms)'})

nexttile(2, [2, 1])
hold on
for idx_loom = 1:20
    bouts_spontaneous = data_parser_new(bouts, 'period', 'bsl', 'window', 'all', 'type', 'immobility', 'nloom', idx_loom);
    [f, x] = ecdf(bouts_spontaneous.durations); 
    plot(x, f, 'Color', col.vars.nloom(10 + idx_loom,:), 'LineWidth', 1);
end

apply_generic(gca, 'xlim', [0 60], 'ylim', [-0.05 1.05], 'font_size', 20)
xlabel('Immobility Frames')
ylabel('ecdf')

exporter(fh, paths_out, 'contaminant_nloom.pdf')