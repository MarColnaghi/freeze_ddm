
clear all

% Load the table first. We will take advantage of an already existing
% dataset.
threshold_imm = 3; threshold_mob = 3; threshold_pc = 4; id_code = sprintf('imm%d_mob%d_pc%d', threshold_imm, threshold_mob, threshold_pc);
paths = path_generator('folder', 'extrema_detection/visual_aid', 'bouts_id', id_code);
load(fullfile(paths.dataset, 'bouts.mat'));
bouts_proc = data_parser_new(bouts, 'type', 'immobility', 'period', 'loom', 'window', 'le');

% Only save useful variables in the table
y = table;
y.sm = bouts_proc.avg_sm_freeze_norm;
y.fs = bouts_proc.avg_fs_1s_norm;
y.ln = bouts_proc.nloom_norm;
y.ls = bouts_proc.sloom_norm;
y.intercept = ones(height(y),1);   % N‑by‑1 column of ones
y.smp = bouts_proc.avg_sm_freeze_norm;
predictors = y.Properties.VariableNames;

points.truncation = [];
points.censoring = 10.5;

fh = figure('color', 'w', 'Position', [100, 100, 800, 300]);
tiledlayout(1, 1, 'TileSpacing', 'compact', 'Padding', 'loose')
nexttile
hold on
apply_generic(gca)

for idx_theta = [0.05 0.1 0.2 0.3]

    bouts = y(1,:);
    params = [idx_theta 0];
    [nll, f, fd] = nll_fly_ddm_newer(params, bouts, points, 'model_simed0', 'iid', 'p', ec);

    plot(fd, f, 'LineWidth', 3, 'DisplayName', sprintf('\\theta: %.2f', idx_theta))
end

xlim([0 1])
ylim([0 0.5])
xlabel('Freeze Duration (s)')
ylabel('pmf')
lh = legend('Box','off', 'FontSize', 18, 'Location', 'westoutside', 'Interpreter','tex');

ax_inset = axes('Position', [0.45 0.5 0.5 0.3]); hold on;

for idx_theta = [0.05 0.1 0.2 0.3]

    bouts = y(1,:);
    params = [idx_theta 0];
    [nll, f, fd] = nll_fly_ddm_newer(params, bouts, points, 'model_simed0', 'iid', 'p', []);

    plot(fd, f, 'LineWidth', 3)
end
xlim([0 11])
ylim([0 0.5])
set(gca, 'YScale', 'log')
apply_generic(gca)
exporter(fh, paths, 'aid_geometric.pdf')