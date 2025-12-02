%% aggregated

id_code = 'imm2_mob2_pc4';
col = cmapper();
paths_out = path_generator('folder', '/spontaneous_process', 'bouts_id', id_code, 'imfirst', false);
bouts = importdata(fullfile(paths_out.dataset, 'bouts.mat'));
bouts_spontaneous = data_parser_new(bouts, 'period', 'bsl', 'window', 'le', 'type', 'immobility', 'nloom', 10:20);
bouts_le = data_parser_new(bouts, 'period', 'loom', 'window', 'le', 'type', 'immobility');

fh = figure('color', 'w', 'Position', [100, 100, 800, 400]);
hold on
histogram(bouts_spontaneous.durations, 'BinMethod', 'integers', 'Normalization', 'probability','EdgeColor','none', 'FaceColor', col.processes.contam)
histogram(bouts_le.durations, 'BinMethod', 'integers', 'Normalization', 'probability', 'EdgeColor','none', 'FaceColor', col.period.loom)
xlims = [-2 302];
ylims = [-0.005 0.205];
apply_generic(gca, 'xlim', xlims, 'ylim', ylims)
xlabel('Duration (frames)')

inset_pos = [0.4, 0.4, 0.2, 0.5];
axes(fh, 'Position', inset_pos)
hold on
apply_generic(gca, 'xlim', xlims, 'ylim', [-0.05 1.05])
ylabel('ecdf')
[f, x] = ecdf(bouts_spontaneous.durations); plot(x, f, 'Color', col.processes.contam, 'LineWidth', 3);
[f, x] = ecdf(bouts_le.durations); plot(x, f, 'Color', col.period.loom, 'LineWidth', 3);

inset_pos(1) = inset_pos(1) + 0.3;
axes(fh, 'Position', inset_pos)
hold on
apply_generic(gca, 'xlim', [-0.03 60.03], 'ylim', [-0.001 0.051])

histogram(bouts_spontaneous.durations, 'BinMethod', 'integers', 'Normalization', 'probability','EdgeColor','none', 'FaceColor', col.processes.contam)
histogram(bouts_le.durations, 'BinMethod', 'integers', 'Normalization', 'probability', 'EdgeColor','none', 'FaceColor', col.period.loom)

% exporter(fh, paths_out, 'bsl_vs_loom_distributions.pdf')


%% Conditions


fh = figure('Position', [100 100 1400 900], 'Color', 'w');
t = tiledlayout(3, 4, 'TileSpacing', 'compact', 'Padding', 'loose');
i = 0;
for idx_sm = 1:3
    for idx_ls = 1:2
        for idx_fs = 1:2

            i = i + 1;
            nexttile(t)
            hold on

            [quant_le, mask] = quantilizer(bouts_le, 'idx_quanti', struct('sm', idx_sm, 'ls', idx_ls, 'fs', idx_fs));
            histogram(quant_le.durations, 1:1:630, 'Normalization', 'probability', 'EdgeColor','none', 'FaceColor', col.period.loom)
            exp_generated = exprnd(1, height(quant_le), 1);
            
            histogram(exp_generated * 60, 1:3:630, 'Normalization', 'probability', 'EdgeColor','none', 'FaceColor', col.period.bsl)

            xlims = [-10 640];
            ylims = [-0.005 0.025];
            apply_generic(gca, 'xlim', xlims, 'ylim', ylims)
            ax = gca;
           % ax.XAxis.Scale = 'log';
        end
    end
end