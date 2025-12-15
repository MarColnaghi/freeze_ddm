clearvars

windows.list = {'cumulative', 'cumulative'};
windows.anchor = {'freeze_onset', 'freeze_onset'};
windows.windows = {'-30', '60'};

col.Set3 = {'#D34E8E', '#31D0CA'};
model_str = 'dddm2';
results_pre = importdata(fullfile('/Users/marcocolnaghi/PhD/freeze_ddm/model_results/causality/fitting_windows/dddm2/cumulative/run01_sizecumul_anchor-loom_onset/-30', sprintf('fit_results_%s.mat', model_str)));
[fh, ax, ax_inset, th] = fd_conditions('results', results_pre, 'no_y', true, 'color', 'gray', 'vis', 'on');
paths = path_generator('folder', 'causality/fitting_windows');

for idx_run = 1:length(windows.list)

    path_code = fullfile(paths.results, ...
        model_str, windows.list{idx_run}, sprintf('run01_sizecumul_anchor-%s', windows.anchor{idx_run}), windows.windows{idx_run});

    results = importdata(fullfile(path_code, sprintf('fit_results_%s.mat', model_str)));
    freezes = importdata(fullfile(path_code, sprintf('surrogate.mat')));

    freezes.sm = freezes.avg_sm_freeze_norm;
    freezes.smp = freezes.avg_sm_pre_norm;
    freezes.fs = freezes.avg_fs_1s_norm;
    freezes.ls = freezes.sloom_norm;
    freezes.ln = freezes.nloom_norm;
    freezes.intercept = ones(height(freezes), 1);
    freezes = freezes(freezes.durations_s >= 0.3, :);


    est_params = table2array(results.estimates_mean(:, ~ismissing(results.estimates_mean)));

    i = 0;
    for idx_sm = 1:3
        for idx_ls = 1:2
            for idx_fs = 1:2

                i = i + 1;
                axes(ax(i));

                [freezes_quant, mask] = quantilizer(freezes, 'idx_quanti', struct('sm', idx_sm, 'ls', idx_ls, 'fs', idx_fs));

                [~, f, fd] =  nll_fly_ddm_newer(est_params, freezes_quant, results.points, strcat('model_', results.fitted_model), 'iid', 'p', []);

                plot(fd, f, 'LineWidth', 1.5, 'Color',  col.Set3{idx_run}, 'LineStyle', '--')
                plot(results.points.censoring, f(end), 'o', 'LineWidth', 1, 'MarkerSize', 5, 'MarkerEdgeColor', col.Set3{idx_run});

                axes(ax_inset(i));
                plot(results.points.censoring, f(end), 'o', 'LineWidth', 1, 'MarkerSize', 5, 'MarkerEdgeColor', col.Set3{idx_run});

            end
        end
    end

end


model_acausal = importdata('/Users/marcocolnaghi/PhD/freeze_ddm/model_results/fitting_freezes/le/dddm2/run11/fit_results_dddm2.mat');
est_params = table2array(model_acausal.estimates_mean(:, ~ismissing(model_acausal.estimates_mean)));

i = 0;
for idx_sm = 1:3
    for idx_ls = 1:2
        for idx_fs = 1:2

            i = i + 1;
            axes(ax(i));

            [freezes_quant, mask] = quantilizer(freezes, 'idx_quanti', struct('sm', idx_sm, 'ls', idx_ls, 'fs', idx_fs));

            [~, f, fd] =  nll_fly_ddm_newer(est_params, freezes_quant, model_acausal.points, strcat('model_', model_acausal.fitted_model), 'iid', 'p', []);

            plot(fd, f, 'k--', 'LineWidth', 1.25)
            axes(ax_inset(i));
            plot(results.points.censoring, f(end), 'ko', 'LineWidth', 1, 'MarkerSize', 5)

        end
    end
end

paths = path_generator('folder', fullfile('causality/fitting_windows', ...
    model_str, windows.list{idx_run}));

export = true;
if export
    exporter(fh, paths, 'fits_xcondition_cmpr.pdf')
end

xl = [0.1 1.5];
yl = [0 3.5];
axes(ax(i));
xlim(xl); ylim(yl);

for idx_tiles = 1:size(th, 1)
    th(idx_tiles, 1).Visible = 'off';
    th(idx_tiles, 2).Visible = 'off';
    th(idx_tiles, 3).Visible = 'off';
    ax(idx_tiles).YAxis.Visible = true;
    axes(ax(idx_tiles))
    apply_generic(gca, 'xticks', xl, 'yticks', yl)
end
exporter(fh, paths, 'fits_cmpr_zoomL.pdf')

xl = [0.5 10.5];
yl = [0 0.5];
xlim(xl); ylim(yl);
for idx_tiles = 1:size(th, 1)
    axes(ax(idx_tiles))
    apply_generic(gca, 'xticks', xl, 'yticks', yl)
end
exporter(fh, paths, 'fits_cmpr_zoomR.pdf')
