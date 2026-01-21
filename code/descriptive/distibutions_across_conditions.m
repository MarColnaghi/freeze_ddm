%% aggregated

id_code = 'imm2_mob2_pc4';
col = cmapper([], 4);
paths_out = path_generator('folder', '/spontaneous_process', 'bouts_id', id_code, 'imfirst', false);
bouts = importdata(fullfile(paths_out.dataset, 'bouts.mat'));

bouts_le = data_parser_new(bouts, 'period', 'loom', 'window', 'le', 'type', 'immobility', 'nloom', 2:20, 'min_dur', 30);

fh = figure('color', 'w', 'Position', [100, 100, 1000, 300]);
tiledlayout(1, 4, 'TileSpacing', 'compact', 'Padding', 'compact')


for idx_sm = 1:4
    
    bout_quant = quantilizer(bouts_le, 'quant', struct('sm', 4, 'fs', 1, 'ls', 2), 'idx_quanti', struct('sm', idx_sm, 'fs', 1, 'ls', 2));

    nexttile
    hold on

    histogram(bout_quant.durations_s, min(bouts_le.durations_s):10/60:max(bouts_le.durations_s), 'Normalization', 'pdf',...
        'EdgeColor', 'none', 'FaceColor', col.vars.sm(idx_sm + 1, :), 'FaceAlpha', 0.35)

    is_censored = bout_quant.durations_s > 10.5;
    find(is_censored)
    [fkde, x] = ksdensity(bout_quant.durations_s, 0:1/60:max(bouts_le.durations_s), ...
        'BoundaryCorrection', 'reflection', 'Bandwidth', 0.16, 'Censoring', is_censored, ...
        'Support', [min(bouts_le.durations_s) - 1e-12, max(bouts_le.durations_s) + 1e-12]);

    plot(x, fkde, 'Color', col.vars.sm(idx_sm + 1, :), 'LineWidth', 1.4)

    apply_generic(gca, 'xlim', [0 10.5], 'ylim', [-0.02 1.52], 'no_y', true)

end