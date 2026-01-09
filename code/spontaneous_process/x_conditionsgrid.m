id_code = 'imm2_mob2_pc4';
col = cmapper();
paths_out = path_generator('folder', '/spontaneous_process', 'bouts_id', id_code, 'imfirst', false);
bouts = importdata(fullfile(paths_out.dataset, 'bouts.mat'));
bouts_spontaneous = data_parser_new(bouts, 'period', 'bsl', 'window', 'le', 'type', 'immobility', 'nloom', 10:20);
bouts_le = data_parser_new(bouts, 'period', 'loom', 'window', 'le', 'type', 'immobility', 'nloom', 2:20);

fh = figure('Position', [100 100 1400 900], 'Color', 'w');
t = tiledlayout(3, 4, 'TileSpacing', 'compact', 'Padding', 'loose');
i = 0;
for idx_sm = 1:3
    for idx_ls = 1:2
        for idx_fs = 1:2

            i = i + 1;
            nexttile(t)
            hold on

            exp_generated = exprnd(1/15, 500000, 1);
            exp_generated = exp_generated(exp_generated >= min(bouts_spontaneous.durations_s));
            histogram(exp_generated * 60, 1:1:630, 'Normalization', 'probability', 'EdgeColor', 'k', 'FaceAlpha', 0.3, 'DisplayStyle', 'stairs')


            [quant_spon, mask] = quantilizer(bouts_spontaneous, 'idx_quanti', struct('sm', idx_sm, 'ls', idx_ls, 'fs', idx_fs));
            histogram(quant_spon.durations, 1:1:630, 'Normalization', 'probability', 'EdgeColor','none', 'FaceColor', col.period.bsl)

%             [quant_le, mask] = quantilizer(bouts_le, 'idx_quanti', struct('sm', idx_sm, 'ls', idx_ls, 'fs', idx_fs));
%             histogram(quant_le.durations, 1:1:630, 'Normalization', 'probability', 'EdgeColor','none', 'FaceColor', col.period.loom)

            xlims = [-10 160];
            ylims = [-0.005 1.025];
            apply_generic(gca, 'xlim', xlims, 'ylim', ylims)
            ax = gca;
            ax.YAxis.Scale = 'log';
        end
    end
end