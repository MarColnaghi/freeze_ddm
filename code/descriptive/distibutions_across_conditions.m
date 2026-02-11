%% aggregated

id_code = 'imm2_mob2_pc4';
paths_out = path_generator('folder', 'descriptive/distributions_across_conditions', 'bouts_id', id_code, 'imfirst', false);
bouts = importdata(fullfile(paths_out.dataset, 'bouts.mat'));

thresholds = define_thresholds;
%thresholds.le_window_fl = [1 60];
%thresholds.le_window_sl = [1 60];
thresholds.le_window_fl = [5 40];
thresholds.le_window_sl = [15 50];

bouts = bouts_formatting(bouts, thresholds);
points.censoring = 10.5;
points.truncation = 0.5;

model_2_fit = 'sddm0';

bouts_le = data_parser_new(bouts, 'period', 'loom', 'window', 'le', 'type', 'immobility', 'nloom', 2:20, 'min_dur', 30);

fh = figure('color', 'w', 'Position', [100, 100, 1200, 400]);
n_quantiles = 4;
col = cmapper([], n_quantiles);

tl = tiledlayout(1, n_quantiles, 'TileSpacing', 'compact', 'Padding', 'compact');
idx_loom = 2;

for idx_sm = 1:n_quantiles

    bout_quant = quantilizer(bouts_le, 'quant', struct('sm', n_quantiles, 'fs', 1), 'idx_quanti', struct('sm', idx_sm, 'fs', 1, 'ls', idx_loom));

    figure(fh)
    tl_h(idx_sm) = nexttile(tl);
    hold on

    histogram(bout_quant.durations_s, min(bouts_le.durations_s):10/60:max(bouts_le.durations_s), 'Normalization', 'pdf',...
        'EdgeColor', 'none', 'FaceColor', col.vars.sm(idx_sm + 1, :), 'FaceAlpha', 0.35)

    is_censored = bout_quant.durations_s > 10.5;
    find(is_censored)
    [fkde, x] = ksdensity(bout_quant.durations_s, 0:1/60:max(bouts_le.durations_s), ...
        'BoundaryCorrection', 'reflection', 'Bandwidth', 0.12, 'Censoring', is_censored, ...
        'Support', [min(bouts_le.durations_s) - 1e-12, max(bouts_le.durations_s) + 1e-12]);

    plot(x, fkde, 'Color', col.vars.sm(idx_sm + 1, :), 'LineWidth', 1.8)


    apply_generic(gca, 'xlim', [0 10.5], 'ylim', [-0.02 1.52], 'no_y', true, 'tick_length', 0.02)
    ax = gca;
    ax.set('Color', 'none')

    model_results = run_fitting_newer(bout_quant, points, model_2_fit, paths_out, 'export', false);
    model_results.quant = quant.(vr);

    est = table2array(model_results.estimates_mean);

    estimates(idx_quantiles, idx_ln, :) = est(~isnan(est));

    bout_quant.sm = bout_quant.avg_sm_freeze_norm;
    bout_quant.smp = bout_quant.avg_sm_freeze_norm;
    bout_quant.tsl = bout_quant.time_since_last;

    bout_quant.fs = bout_quant.avg_fs_1s_norm;
    bout_quant.ls = bout_quant.sloom_norm;
    bout_quant.ln = bout_quant.nloom_norm;
    bout_quant.intercept = ones(height(bout_quant), 1);

    [~, f, fd] = nll_fly_ddm_newer( est(~isnan(est)), bout_quant, points, strcat('model_', model_2_fit), 'iid', 'p');

    plh(idx_sm) = plot(fd, f, 'LineWidth', 1.8, 'LineStyle', '--', 'Color', 'k');
end

%%
xlabel('Freeze Duration (s)')
exporter(fh, paths_out, sprintf('across_sm_%d.pdf', idx_loom))

figure(fh)

for idx_sm = 1:n_quantiles
    
    axes(tl_h(idx_sm))

end