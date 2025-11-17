
model = 'simed0';
run = 'run03';
paths = path_generator('folder', fullfile('fitting_freezes/le/quantiles', model, run));

var = 'avg_sm_freeze_norm';
%var = 'avg_fs_1s_norm';
%var = 'nloom_norm';

estimates = importdata(fullfile(paths.results, var, 'estimates.mat'));
model_results = importdata(fullfile(paths.results, var, 'model_results.mat'));

fh = figure('color','w','Position',[100, 100, 700, 300]);
tl = tiledlayout(1, size(estimates, 3), 'TileSpacing', 'compact', 'Padding', 'compact')

x = 1:size(estimates, 1);

for idx_params = 1:size(estimates, 3)
    nexttile
    hold on
    scatter(x, estimates(:, 1, idx_params), 75, 'r', 'filled')
    scatter(x, estimates(:, 2, idx_params), 75, 'b', 'filled')
    axis square

    xlim([0.5 4.5])
    xticks(x)
    xticklabels({'Q1', 'Q2', 'Q3', 'Q4'})

    if idx_params == 1
        ylim([0.15 0.4])

    elseif idx_params == 2
        ylim([0 0.5])
    end

    if strcmp(var, 'avg_sm_freeze_norm')
        xlabel('Social Motion Range');
    elseif strcmp(var, 'avg_fs_1s_norm')
        xlabel('Focal Speed Range');
    elseif strcmp(var, 'nloom_norm')
        xlabel('Loom Number');
    end

    apply_generic(gca)
end

exporter(fh, paths, sprintf('%s.pdf', var))