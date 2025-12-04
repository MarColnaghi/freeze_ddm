% Double DDM test

clearvars

% Load the table first. We will take advantage of an already existing
% dataset.
threshold_imm = 2; threshold_mob = 2; threshold_pc = 4; id_code = sprintf('imm%d_mob%d_pc%d', threshold_imm, threshold_mob, threshold_pc);
model_str = 'dddm2';
col.Set3 = cbrewer2('Set3', 5);
    

fh = figure('Position', [100 100 900 600], 'Color', 'w');
t = tiledlayout(1, 1, 'TileSpacing', 'loose', 'Padding', 'loose');
mkrsize = 220;
nexttile

for idx_run = 1:2
    run_str = sprintf('run0%d', idx_run);

    str_folder = fullfile('causality/fitting_windows', model_str, run_str);
    paths = path_generator('folder', str_folder, 'bouts_id', id_code, 'imfirst', false);
    d = dir(paths.results);
    d = d(~ismember({d.name}, {'.','..', '.DS_Store'}));

    elbos = nan(size(d, 1), 2);
    get_window_size = nan(size(d, 1), 1);

    for idx_windows = 1:size(d, 1)
        get_window_size(idx_windows) = str2double(d(idx_windows).name);
        model_results = importdata(fullfile(d(idx_windows).folder, d(idx_windows).name, sprintf('fit_results_%s.mat', model_str)));
        elbos(idx_windows, :) = model_results.elbo;
        
    end

    hold on

    scatter(-get_window_size, elbos(:,1), mkrsize, 'o', 'MarkerFaceColor', col.Set3(2 + idx_run, :), ...
        'MarkerEdgeColor', 'none')
    [sorted_window_size, i] = sort(-get_window_size);
    plot(sorted_window_size, elbos(i,1), 'Color', col.Set3(2 + idx_run, :), 'LineStyle', '--')
end

model_acausal = importdata('/Users/marcocolnaghi/PhD/freeze_ddm/model_results/fitting_freezes/le/dddm2/run09/fit_results_dddm2.mat');
elbo_acausal = model_acausal.elbo;

scatter(30, elbo_acausal(:,1), mkrsize, 'o', 'MarkerFaceColor', 'k', ...
        'MarkerEdgeColor', 'none')
apply_generic(gca, 'xlim', [-610 60], 'xticks', sort([-get_window_size(:); 30]))
xlabel('Integration Windows (frames)')
ylabel('ELBO')
xline(0, 'k-', 'LineWidth',3)

exporter(fh, paths, 'elbo_fitted_models.pdf')