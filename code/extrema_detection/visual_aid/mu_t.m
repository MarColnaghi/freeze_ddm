
clear all

% Load the table first. We will take advantage of an already existing
% dataset.
threshold_imm = 3; threshold_mob = 3; threshold_pc = 4; id_code = sprintf('imm%d_mob%d_pc%d', threshold_imm, threshold_mob, threshold_pc);
paths = path_generator('folder', 'extrema_detection/visual_aid', 'bouts_id', id_code);
load(fullfile(paths.dataset, 'bouts.mat'));
bouts_proc = data_parser_new(bouts, 'type', 'immobility', 'period', 'loom', 'window', 'le');
motion_cache = importdata(fullfile(paths.cache_path, 'motion_cache.mat'));

%  Now we added our vector column to the bouts table

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

chunk_len = points.censoring * 60;

for idx_trials = 1:height(bouts_proc)

    ons = bouts_proc.onsets(idx_trials);
    sum_motion = motion_cache(bouts_proc.fly(idx_trials));
    sm_raw{idx_trials} = sum_motion(ons:ons + chunk_len - 1) ./ 10;

end

soc_mot_array = cell2mat(sm_raw)';
extra.soc_mot_array = soc_mot_array;

%%
col.order = cbrewer2('Set2', 8);

fh = figure('color', 'w', 'Position', [100, 100, 700, 300]);
tiledlayout(1, 1, 'TileSpacing', 'compact', 'Padding', 'loose')
nexttile
hold on
ax = gca;
apply_generic(gca)
idx_trial = 19;
i = 0;
theta = 0.3;

for idx_mu = [15, 10, 5, 0]

    i = i + 1;
    bouts = y(idx_trial,:);
    params = [idx_mu theta 0];
    ec.soc_mot_array = extra.soc_mot_array(idx_trial, :);
    [nll, f, fd] = nll_fly_ddm_newer(params, bouts, points, 'model_ed1', 'iid', 'p', ec);

    plot(fd, f + 0.25 * i, 'Marker', 'o', 'LineWidth', 2, 'DisplayName', sprintf('\\beta_{sm}: %.1f', idx_mu), 'Color', col.order(i, :))
    plot([-0.05 -0.05], [0.25 * i; 0.25 * i + 0.1], 'LineWidth', 2, 'Color', 'k', 'HandleVisibility', 'off')
end
plot(0, 0, 'Visible', 'off', 'DisplayName', sprintf('\\theta: %.1f', theta), 'Color', 'w')

text(-0.065, (0.25 * i + 0.25 * i + 0.1) / 2, '0.1', 'VerticalAlignment', 'middle', 'HorizontalAlignment','right', 'FontSize', 24)
imagesc(fd, -0.3, ec.soc_mot_array)
ax.YAxis.Visible = 'off';
xticks([]);
colormap(cbrewer2('Reds'))
ylim([0.1 1.3])
xlim([-0.1 1.2])
xlabel('samples')
ylabel('pmf')
lh = legend('Box','off', 'FontSize', 18, 'Location', 'westoutside', 'Interpreter','tex');
apply_generic(gca)
colororder();

exporter(fh, paths, 'mu_sensitivity.pdf')
