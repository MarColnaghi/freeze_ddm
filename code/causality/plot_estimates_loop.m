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

for idx_run = 1
    run_str = sprintf('run0%d', idx_run);

    str_folder = fullfile('causality/fitting_windows', model_str, run_str);
    paths = path_generator('folder', str_folder, 'bouts_id', id_code, 'imfirst', false);
    d = dir(paths.results);
    d = d(~ismember({d.name}, {'.','..', '.DS_Store'}));


    for idx_windows = 1:size(d, 1)
        model_results = importdata(fullfile(d(idx_windows).folder, d(idx_windows).name, sprintf('fit_results_%s.mat', model_str)));
        fh = plot_estimates('results', model_results, 'export', false, 'ylimits', [-2 4]);
        paths_window.fig = fullfile(paths.fig, d(idx_windows).name);
        exporter(fh, paths_window, 'estimates.pdf')
    end
end

%%
fh = figure('Position', [100 100 800 310], 'Color', 'w');
t = tiledlayout(1, 1, 'TileSpacing', 'loose', 'Padding', 'loose');
col.Set3 = cbrewer2('Set3', 9);
    hold on

for idx_run = 1
    run_str = sprintf('run0%d', idx_run);

    str_folder = fullfile('causality/fitting_windows', model_str, run_str);
    paths = path_generator('folder', str_folder, 'bouts_id', id_code, 'imfirst', false);
    d = dir(paths.results);
    d = d(~ismember({d.name}, {'.','..', '.DS_Store'}));

    estimates = nan(10, 1);
    get_window_size = nan(size(d, 1), 1);

    for idx_windows = 1:size(d, 1)
        get_window_size(idx_windows) = str2double(d(idx_windows).name);
        model_results = importdata(fullfile(d(idx_windows).folder, d(idx_windows).name, sprintf('fit_results_%s.mat', model_str)));
        est_temp = table2array(model_results.estimates_mean);
        est_temp = est_temp(~isnan(est_temp));
        estimates(:) = est_temp;
        plot(1:10, estimates, 'o', 'Color', col.Set3(idx_windows,:), 'MarkerFaceColor', col.Set3(idx_windows + 2,:), ...
            'MarkerSize', 7, 'DisplayName', d(idx_windows).name, 'MarkerEdgeColor', 'k');
    end

end

model_acausal = importdata('/Users/marcocolnaghi/PhD/freeze_ddm/model_results/fitting_freezes/le/dddm2/run09/fit_results_dddm2.mat');
est_temp = table2array(model_acausal.estimates_mean);
est_temp = est_temp(~isnan(est_temp));

plot(1:10, estimates, 'ko', 'LineWidth', 1.2, 'MarkerFaceColor', 'k', ...
            'MarkerSize', 8, 'DisplayName', 'acausal');
lh = legend('Location', 'eastoutside', 'Box', 'off', 'FontSize', 22);
apply_generic(gca, 'xticks', 1:10);
hold on

est_means = model_results.estimates_mean(:, ~ismissing(model_results.estimates_mean));
est_std = model_results.estimates_std(:, ~ismissing(model_results.estimates_mean));
[suffixes, prefixes] = extract_dep(est_means);

% Replace 'intercept' with '0' in the suffixes
suffixes_replaced = suffixes;
suffixes_replaced(strcmp(suffixes, 'intercept')) = {'0'};

% Combine using cellfun
xx = 1:size(suffixes, 2);
ax = gca;
result = cellfun(@(pre, suf) ['$$\beta_{' pre '}^{' suf '}$$'], ...
                 prefixes, suffixes_replaced, 'UniformOutput', false);
xticklabels(result); set(ax.XAxis, 'TickLabelInterpreter', 'latex', 'FontSize', 24);


exporter(fh, paths, 'estimates_change.pdf')
function [suffixes, prefixes] = extract_dep(results)

params = results.Properties.VariableNames;
parts_str = cellfun(@(s) split(s, '_'), params, 'UniformOutput', false);
suffixes = cellfun(@(parts) parts{end}, parts_str, 'UniformOutput', false);
prefixes = cellfun(@(concat) concat{1:end-1}, parts_str, 'UniformOutput', false);
end
