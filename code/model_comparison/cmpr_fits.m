clearvars

% Load the table first. We will take advantage of an already existing
% dataset.
threshold_imm = 2; threshold_mob = 2; threshold_pc = 4; id_code = sprintf('imm%d_mob%d_pc%d', threshold_imm, threshold_mob, threshold_pc);
model_str = 'dddm2';
model_func = str2func(strcat('model_',model_str));
model_struct = model_func();

col = cmapper();

color_list = {col.timevarying_sm, col.extremadetection};
model_list = {'dddm2', 'ded2'};
run_list = {'run13', 'run10'};
label_list = {'Integration', 'Extrema Detection'};

paths = path_generator('folder', 'fitting_freezes/le');

results_folder = fullfile(paths.results, model_list{1}, run_list{1});
temp_results = importdata(fullfile(results_folder, sprintf('fit_results_%s.mat', model_list{1})));
surrogate = importdata(fullfile(results_folder, 'surrogate.mat'));

[fh, ax, ax_inset, th] = fd_conditions('results', temp_results, 'no_y', true, 'col', 'gray', 'vis', 'on');

for idx_model = 1:length(model_list)

    results_folder = fullfile(paths.results, model_list{idx_model}, run_list{idx_model});

    model_results = importdata(fullfile(results_folder, sprintf('fit_results_%s.mat', model_list{idx_model})));
    freezes = importdata(fullfile(results_folder, 'surrogate.mat'));
    extra = importdata(fullfile(results_folder, 'extra.mat'));

    est_params = table2array(model_results.estimates_mean(:, ~ismissing(model_results.estimates_mean)));

    model_results.points.truncation

    i = 0;
    for idx_sm = 1:3
        for idx_ls = 1:2
            for idx_fs = 1:2

                i = i + 1;
                axes(ax(i));


                [freezes_quant, mask] = quantilizer(freezes, 'idx_quanti', struct('sm', idx_sm, 'ls', idx_ls, 'fs', idx_fs));

                ec = extra;
                ec.soc_mot_array = extra.soc_mot_array(mask, :);
                ec = rmfield(ec, 'tndt');
                [~, f, fd] =  nll_fly_ddm_newer(est_params, freezes_quant, model_results.points, strcat('model_', model_results.fitted_model), 'iid', 'p', ec);

                plot(fd, f, 'LineWidth', 1.5, 'Color',  color_list{idx_model}, 'LineStyle', '--')
                plot(model_results.points.censoring, f(end), 'o', 'LineWidth', 1, 'MarkerSize', 5, 'MarkerEdgeColor', color_list{idx_model});

                axes(ax_inset(i));
                plot(model_results.points.censoring, f(end), 'o', 'LineWidth', 1, 'MarkerSize', 5, 'MarkerEdgeColor', color_list{idx_model});

            end
        end
    end
end

paths = path_generator('folder', 'model_comparison');
exporter(fh, paths, 'cmpr_fits.pdf')

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