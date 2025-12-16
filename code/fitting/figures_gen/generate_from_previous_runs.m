export = false;
paths = path_generator('folder', 'fitting_freezes/le');
model = 'ded2';
run = 'run10';
results_folder = fullfile(paths.results, model, run);

model_results = importdata(fullfile(results_folder, sprintf('fit_results_%s.mat', model)));
freezes = importdata(fullfile(results_folder, 'surrogate.mat'));
% plot_estimates('results', model_results, 'export', export, 'ylimits', [-2 5])
[fh, ax, ax_inset] = fd_conditions('results', model_results, 'no_y', true);
overlay_fits(fh, ax, ax_inset, 'results', model_results, 'export', true, 'extra', extra)
