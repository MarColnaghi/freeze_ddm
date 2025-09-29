% This code shows how a stationary drift ddm works.

close all
clearvars
clc
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
plot(x_plot(60:end-60), sm_aroundfreeze(60:end-60), 'Color', 'w', 'LineWidth', 3 )
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

%export
%%
plot([mean(sm_aroundfreeze(60:end-60)) mean(sm_aroundfreeze(60:end-60))], [0 mean(sm_aroundfreeze(60:end-60)) * gt(idx)], 'Color', col.stationary_sm, 'LineStyle', ':', 'LineWidth', 4)
plot([0 mean(sm_aroundfreeze(60:end-60))], [mean(sm_aroundfreeze(60:end-60)) * gt(idx) mean(sm_aroundfreeze(60:end-60)) * gt(idx)], 'Color', col.stationary_sm, 'LineStyle', ':', 'LineWidth', 4)
sch = scatter(0, mean(sm_aroundfreeze(60:end-60)) * gt(idx), 500, 'filled', 'Clipping', 'off', 'MarkerFaceColor', col.visual_aid.mu);
set(gca ,'Layer', 'Top')


%%
sim_params.dt = 1/60;
sim_params.T = 30;
sim_params.time_vector = 0:sim_params.dt:sim_params.T;
sim_params.z = 0;

fh = figure('color','w', 'Position', [100 100 700 400]);
tiledlayout(2, 1, 'TileSpacing', 'compact', 'Padding', 'compact')

nexttile(1, [1,1])
hold on
model_eval = eval(sprintf('model_%s', model));
ncomp_vars = evaluate_model(model_eval, gt, bout_tbl);
rt_st = nan(1,1);
for idx_trials = 1
    [rt_st(idx_trials), traj_st] = drift_diff_new('mu', ncomp_vars.mu2, 'theta', 0, ...
        'z', -ncomp_vars.theta2, 'dt', sim_params.dt, 'T', sim_params.T, 'ndt',  ncomp_vars.tndt);
    plot(traj_st, 'k-', 'LineWidth', 1)
end
xlim([0 30*60])
yticks([])
rt_st = rt_st * 60;
ax2 = gca;
ax2.XAxisLocation = 'origin';
 xticks([])
ylabel('Safety Evidence')
xline(0,  'k--', 'LineWidth', 2)

apply_generic(ax2, 18);
% nexttile(1)
% histogram(rt_st, 50, 'Normalization', 'pdf')
% xlim([0 30*60])
% ax = gca;
% yticks([])
% apply_generic(ax)

nexttile(2)
hold on
x_plot = (bout_tbl.onsets - 60: bout_tbl.onsets + bout_tbl.durations + 60) - bout_tbl.onsets;
plot(x_plot, sm_aroundfreeze, 'k-', 'LineWidth', 2)
plot(x_plot(60:end-60), sm_aroundfreeze(60:end-60), 'Color', 'w', 'LineWidth', 3 )
plot(x_plot(60:end-60), repmat(mean(sm_aroundfreeze(60:end-60)), length(x_plot(60:end-60)), 1), 'Color', col.stationary_sm, 'LineWidth', 3)
ax3 = gca;
ax3.XAxis.Visible = 'off';
apply_generic(ax3)
ylabel('Social Motion')
xline(0,  'k--', 'Label', 'Freeze Entry', 'LineWidth', 2, 'FontSize', 20, 'LabelOrientation', 'horizontal')
xline(bout_tbl.durations, 'k--', 'Label', 'Freeze Exit', 'LineWidth', 2, 'FontSize', 20, 'LabelOrientation', 'horizontal')
ylim([0 4])
yticks([])
apply_generic(ax3, 18);
linkaxes([ax2, ax3], 'x')
xlim([-60 660])


%%
nexttile(1, [1,1])
hold on

model_eval = eval(sprintf('model_%s', model));
ncomp_vars = evaluate_model(model_eval, gt, bout_tbl);

% Run a single simulated trajectory (as you had)
[rt_sample, traj_st] = drift_diff_new('mu', ncomp_vars.mu2, 'theta', 0, ...
    'z', -ncomp_vars.theta2, 'dt', sim_params.dt, 'T', sim_params.T, 'ndt', ncomp_vars.tndt);

% Pre-create graphics objects for efficiency
hLine = plot(nan, nan, 'k-', 'LineWidth', 1);                      % growing trajectory
hDot  = plot(1, traj_st(1), 'ko', 'MarkerFaceColor','k', 'MarkerSize', 6);  % moving dot

% Choose an arrow length relative to the y-range; update after initial limits set
% (If you want a fixed length, set arrowLen to a constant.)
arrowLen = 2;  % temporary; will rescale after ylim known
hArrow = quiver(2, traj_st(1), 3, arrowLen, 0, ...   % last arg 0 = no autoscale
    'MaxHeadSize', 2, 'Color', 'k', 'LineWidth', 1.5);

% Axis cosmetics (as you had)
yticks([])
ax2 = gca;
ax2.XAxisLocation = 'origin';
xticks([])
ylabel('Safety Evidence')
xline(0,  'k--', 'LineWidth', 2)

apply_generic(ax2, 18);

% Make sure the y-limits are set before choosing arrow length
drawnow;
yl = ylim;
arrowLen = 1; 
set(hArrow, 'UData', ncomp_vars.mu2, 'VData', arrowLen);  % point arrow upward

% --- Animation loop ---
n = numel(traj_st(~isnan(traj_st)));
% Optional: animate ~real-time; set to sim_params.dt for 1x speed, or smaller for faster
framePause = sim_params.dt;   % change to, e.g., 0 for fastest rendering

for k = 2:n
    % Update line to current time
    set(hLine, 'XData', 1:k, 'YData', traj_st(1:k));

    % Move dot
    set(hDot, 'XData', k, 'YData', traj_st(k));

    % Keep arrow at the dot, pointing up
    set(hArrow, 'XData', k, 'YData', traj_st(k), 'UData', 5, 'VData', 1.2);

    % Render
    drawnow limitrate
    if framePause > 0
        pause(framePause);
    end
end

% Convert RT to seconds if you still need it later
rt_st = rt_sample * 60;
