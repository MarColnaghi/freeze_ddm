
col.Set3 = cbrewer2('Set3');
model_str = 'dddm2';
results = importdata(fullfile('/Users/marcocolnaghi/PhD/freeze_ddm/model_results/causality/fitting_windows/dddm2/run01/60', sprintf('fit_results_dddm2.mat', model_str)));
[fh, ax, ax_inset] = fd_conditions('results', results, 'no_y', true, 'color', 'gray');
selected_windows = {'60'};

for idx_run = 1:2
    run_str = sprintf('run0%d', idx_run);

    str_folder = fullfile('causality/fitting_windows', model_str, run_str);
    paths = path_generator('folder', str_folder, 'imfirst', false);
    d = dir(paths.results);
    d = d(~ismember({d.name}, {'.','..', '.DS_Store'}));windows = [600, 300, 60];

    freezes = importdata(fullfile(d(1).folder, d(1).name, 'surrogate.mat'));

    freezes.sm = freezes.avg_sm_freeze_norm;
    freezes.smp = freezes.avg_sm_pre_norm;
    freezes.fs = freezes.avg_fs_1s_norm;
    freezes.ls = freezes.sloom_norm;
    freezes.ln = freezes.nloom_norm;
    freezes.intercept = ones(height(freezes), 1);
    freezes = freezes(freezes.durations_s >= 0.3, :);

    for idx_window = selected_windows
        window = idx_window{1};
        idx_folder_list = find(strcmp(window, {d.name}));
        results = importdata(fullfile(d(idx_folder_list).folder, d(idx_folder_list).name, sprintf('fit_results_%s.mat', model_str)));
        est_params = table2array(results.estimates_mean(:, ~ismissing(results.estimates_mean)));

        i = 0;
        for idx_sm = 1:3
            for idx_ls = 1:2
                for idx_fs = 1:2

                    i = i + 1;
                    axes(ax(i));

                    [freezes_quant, mask] = quantilizer(freezes, 'idx_quanti', struct('sm', idx_sm, 'ls', idx_ls, 'fs', idx_fs));

                    [~, f, fd] =  nll_fly_ddm_newer(est_params, freezes_quant, results.points, strcat('model_', results.fitted_model), 'iid', 'p', []);

                    plot(fd, f, 'LineWidth', 1.25, 'Color', col.Set3(2 + idx_run, :), 'LineStyle', '--')

                    axes(ax_inset(i));
                    plot(results.points.censoring, f(end), 'o', 'LineWidth', 1, 'MarkerSize', 5, 'MarkerEdgeColor', col.Set3(2 + idx_run,:));

                end
            end
        end

    end


end

model_acausal = importdata('/Users/marcocolnaghi/PhD/freeze_ddm/model_results/fitting_freezes/le/dddm2/run09/fit_results_dddm2.mat');
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


if export
    paths.fig = results.fig_path;
    if conditions
        exporter(fh, paths, 'fits_xcondition.pdf')

    else
        exporter(fh, paths, 'fits.pdf')
    end
end
