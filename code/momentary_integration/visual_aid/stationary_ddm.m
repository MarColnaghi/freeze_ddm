% This code shows how a stationary drift ddm works.

threshold_imm = 3; threshold_mob = 3; threshold_pc = 4; id_code = sprintf('imm%d_mob%d_pc%d', threshold_imm, threshold_mob, threshold_pc);

exporting = true;
col = cmapper();
sim_params.rng = 1;
rng(sim_params.rng);

% Model
model = 'dddm2';
select_run = 'run03';

% Load important files
paths = path_generator('folder', fullfile('/fitting_freezes/le', model, select_run));
load(fullfile(paths.results, 'surrogate.mat'))
load(fullfile(paths.results, sprintf('fit_results_%s.mat', model)));
gt = table2array(estimates_mean(:, find(~ismissing(estimates_mean))));
gt_table = estimates_mean(:, find(~ismissing(estimates_mean)));
surrogate = surrogate;

% Load the motion ts
sim_params.motion_cache_path = fullfile(paths.dataset, 'motion_cache.mat');
load(sim_params.motion_cache_path)
paths = path_generator('folder', 'momentary_evidence/visual_aid', 'bouts_id', id_code);

% Plot Drift Coefficient
y = table;
y.sm = surrogate.avg_sm_freeze_norm;
y.fs = surrogate.avg_fs_1s_norm;
y.ln = surrogate.nloom_norm;
y.ls = surrogate.sloom_norm;
y.intercept = ones(height(y),1);
y.smp = surrogate.avg_sm_freeze_norm;
y.onsets = surrogate.onsets;
y.fly = surrogate.fly;
y.id = surrogate.id;
y.durations = surrogate.durations;

%%

for idx_bout = 1120
bout_tbl = y(idx_bout, :);
sm = motion_cache(bout_tbl.fly);
sm_aroundfreeze = sm(bout_tbl.onsets - 60: bout_tbl.onsets + bout_tbl.durations + 60);

fh = figure('color','w', 'Position', [100 100 1000 250]);
hold on
x_plot = (bout_tbl.onsets - 60: bout_tbl.onsets + bout_tbl.durations + 60) - bout_tbl.onsets;
plot(x_plot, sm_aroundfreeze, 'k-', 'LineWidth', 2)
plot(x_plot(60:end-60), sm_aroundfreeze(60:end-60), 'Color', col.timevarying_sm, 'LineWidth', 3)
ax = gca;
ax.XAxis.Visible = 'off';
apply_generic(ax)
ylabel('Social Motion')
xline(0,  'k--', 'Label', 'Freeze Entry', 'LineWidth', 2, 'FontSize', 20, 'LabelOrientation', 'horizontal')
xline(bout_tbl.durations, 'k--', 'Label', 'Freeze Exit', 'LineWidth', 2, 'FontSize', 20, 'LabelOrientation', 'horizontal')
ylim([0 4])
end

for idx_bout = 1120
bout_tbl = y(idx_bout, :);
sm = motion_cache(bout_tbl.fly);
sm_aroundfreeze = sm(bout_tbl.onsets - 60: bout_tbl.onsets + bout_tbl.durations + 60);

fh = figure('color','w', 'Position', [100 100 1000 250]);
hold on
x_plot = (bout_tbl.onsets - 60: bout_tbl.onsets + bout_tbl.durations + 60) - bout_tbl.onsets;
plot(x_plot, sm_aroundfreeze, 'k-', 'LineWidth', 2)
plot(x_plot(60:end-60), sm_aroundfreeze(60:end-60), 'Color', [hex2rgb(col.timevarying_sm), 0.4], 'LineWidth', 3 )
plot(x_plot(60:end-60), repmat(mean(sm_aroundfreeze(60:end-60)), length(x_plot(60:end-60)), 1), 'Color', col.stationary_sm, 'LineWidth', 3)
ax = gca;
ax.XAxis.Visible = 'off';
apply_generic(ax)
ylabel('Social Motion')
xline(0,  'k--', 'Label', 'Freeze Entry', 'LineWidth', 2, 'FontSize', 20, 'LabelOrientation', 'horizontal')
xline(bout_tbl.durations, 'k--', 'Label', 'Freeze Exit', 'LineWidth', 2, 'FontSize', 20, 'LabelOrientation', 'horizontal')
ylim([0 4])


end

%%

idx = find(contains(string(gt_table.Properties.VariableNames), 'mu1'));     % -> likely returns []

xxtick = [0 , 1];

fh = figure('color','w', 'Position', [100 100 500 400]);
hold on

plot(xxtick, xxtick * gt(idx), 'LineWidth', 10, 'Color', col.vars.sm(130,:))
ax = gca;
set(ax,'linewidth', 7,'TickDir','both');
set(ax, 'box','off')
yticks(0);
xticks([]);
xlim([0, 0.5])
ylim([-0.5, 0.5])

ax.XAxisLocation = 'origin';
%set(gca ,'Layer', 'Top')

xlabel('Social Motion')
ylabel('Drift Rate (\mu)', 'Interpreter', 'tex', 'Rotation', 90)
ax.XAxis.FontSize = 28;
ax.YAxis.FontSize = 32;
set(ax, 'color', 'None');
set(ax,'TickLength',[0.02, 0.02])



plot([mean(sm_aroundfreeze(60:end-60)) mean(sm_aroundfreeze(60:end-60))], [0 mean(sm_aroundfreeze(60:end-60)) * gt(idx)], 'Color', col.stationary_sm, 'LineStyle', ':', 'LineWidth', 4)
plot([0 mean(sm_aroundfreeze(60:end-60))], [mean(sm_aroundfreeze(60:end-60)) * gt(idx) mean(sm_aroundfreeze(60:end-60)) * gt(idx)], 'Color', col.stationary_sm, 'LineStyle', ':', 'LineWidth', 4)
sch = scatter(0, mean(sm_aroundfreeze(60:end-60)) * gt(idx), 500, 'filled', 'Clipping', 'off', 'MarkerFaceColor', col.visual_aid.mu);
set(gca ,'Layer', 'Top')


%%
sim_params.dt = 1/60;
sim_params.T = 30;
sim_params.time_vector = 0:sim_params.dt:sim_params.T;
sim_params.z = 0;

fh = figure('color','w', 'Position', [100 100 1000 900]);
tiledlayout(3, 1, 'TileSpacing', 'tight', 'Padding', 'tight')

nexttile(2, [1,1])
hold on
model_eval = eval(sprintf('model_%s', model));
ncomp_vars = evaluate_model(model_eval, gt, bout_tbl);
rt_st = nan(1,1);
for idx_trials = 1
    [rt_st(idx_trials), traj_st] = drift_diff_new('mu', ncomp_vars.mu2, 'theta', 0, ...
        'z', -ncomp_vars.theta2, 'dt', sim_params.dt, 'T', sim_params.T, 'ndt',  ncomp_vars.tndt);
    plot(traj_st, 'k-', 'LineWidth', 0.2)
end
xlim([0 30*60])
yticks([])
rt_st = rt_st * 60;
ax = gca;
ax.XAxisLocation = 'origin';
apply_generic(ax); xticks([])

nexttile(1)
histogram(rt_st, 50, 'Normalization', 'pdf')
xlim([0 30*60])
ax = gca;
yticks([])
apply_generic(ax)

nexttile(3)
hold on
x_plot = (bout_tbl.onsets - 60: bout_tbl.onsets + bout_tbl.durations + 60) - bout_tbl.onsets;
plot(x_plot, sm_aroundfreeze, 'k-', 'LineWidth', 2)
plot(x_plot(60:end-60), sm_aroundfreeze(60:end-60), 'Color', [hex2rgb(col.timevarying_sm), 0.4], 'LineWidth', 3 )
plot(x_plot(60:end-60), repmat(mean(sm_aroundfreeze(60:end-60)), length(x_plot(60:end-60)), 1), 'Color', col.stationary_sm, 'LineWidth', 3)
ax = gca;
ax.XAxis.Visible = 'off';
apply_generic(ax)
ylabel('Social Motion')
xline(0,  'k--', 'Label', 'Freeze Entry', 'LineWidth', 2, 'FontSize', 20, 'LabelOrientation', 'horizontal')
xline(bout_tbl.durations, 'k--', 'Label', 'Freeze Exit', 'LineWidth', 2, 'FontSize', 20, 'LabelOrientation', 'horizontal')
ylim([0 4])

%%
axes
x_plot = (bout_tbl.onsets - 60: bout_tbl.onsets + bout_tbl.durations + 60) - bout_tbl.onsets;
plot(x_plot, sm_aroundfreeze, 'k-', 'LineWidth', 2)
plot(x_plot(60:end-60), sm_aroundfreeze(60:end-60), 'Color', col.timevarying_sm, 'LineWidth', 3)
ax = gca;
ax.XAxis.Visible = 'off';
apply_generic(ax)
ylabel('Social Motion')
xline(0,  'k--', 'Label', 'Freeze Entry', 'LineWidth', 2, 'FontSize', 20, 'LabelOrientation', 'horizontal')
xline(bout_tbl.durations, 'k--', 'Label', 'Freeze Exit', 'LineWidth', 2, 'FontSize', 20, 'LabelOrientation', 'horizontal')
ylim([0 4])