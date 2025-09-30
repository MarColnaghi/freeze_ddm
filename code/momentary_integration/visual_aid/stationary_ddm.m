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

% Load the motion ts
sim_params.motion_cache_path = fullfile(paths.dataset, 'motion_cache.mat');
load(sim_params.motion_cache_path)
paths = path_generator('folder', 'momentary_integration/visual_aid', 'bouts_id', id_code);

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

pick_bout = 1120;
pad = 60;

for idx_bout = pick_bout
    bout_tbl = y(idx_bout, :);
    ons = bout_tbl.onsets;
    dur = bout_tbl.durations;
    sm = motion_cache(bout_tbl.fly);

    sm_aroundfreeze = sm(ons - pad: ons + dur + pad);

    fh = figure('color','w', 'Position', [100 100 1000 250]);
    tl = tiledlayout(1, 1, 'TileSpacing', 'compact', 'Padding', 'loose');

    nexttile
    hold on
    x_plot = (ons - pad: ons + dur + pad) - ons;
    plot(x_plot, sm_aroundfreeze, 'k-', 'LineWidth', 2)
    lh = plot(x_plot(pad:end-pad), sm_aroundfreeze(pad:end-pad), 'Color', col.timevarying_sm, 'LineWidth', 3);
    ax = gca;
    ax.XAxis.Visible = 'off';

    apply_generic(ax)
    ylabel('Social Motion')
    xline(0,  'k--', 'Label', 'Freeze Entry', 'LineWidth', 2, 'FontSize', 20, 'LabelOrientation', 'horizontal')
    xline(bout_tbl.durations, 'k--', 'Label', 'Freeze Exit', 'LineWidth', 2, 'FontSize', 20, 'LabelOrientation', 'horizontal')
    %yticks([ax.YTick 5.44])
    ylim([0 4])

    exporter(fh, paths, 'tv_signal.pdf', 'flag', true )

    lh.Color = 'w';
    figure(fh)
    plot(x_plot(pad:end-pad), repmat(mean(sm_aroundfreeze(pad:end-pad)), length(x_plot(pad:end-pad)), 1), 'Color', col.stationary_sm, 'LineWidth', 3)
    %yticks(sort([ax.YTick round(mean(sm_aroundfreeze(pad:end-pad)),2)]))
    exporter(fh, paths, 'st_signal.pdf', 'flag', true )

end


%% Now let's plot the relationship

idx = find(contains(string(gt_table.Properties.VariableNames), 'mu1')); 

xxtick = [0 , 1];

fh = figure('color','w', 'Position', [100 100 500 400]);
tl = tiledlayout(1, 1, 'TileSpacing', 'compact', 'Padding', 'loose');
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

xlabel('Social Motion')
ylabel('Drift Rate (\mu)', 'Interpreter', 'tex', 'Rotation', 90)
ax.XAxis.FontSize = 28;
ax.YAxis.FontSize = 32;
set(ax, 'color', 'None');
set(ax,'TickLength',[0.02, 0.02])

exporter(fh, paths, 'slope_sm.pdf', 'flag', true )

plot([mean(sm_aroundfreeze(pad:end-pad)) mean(sm_aroundfreeze(pad:end-pad))], [0 mean(sm_aroundfreeze(pad:end-pad)) * gt(idx)], 'Color', col.stationary_sm, 'LineStyle', ':', 'LineWidth', 4)
plot([0 mean(sm_aroundfreeze(pad:end-pad))], [mean(sm_aroundfreeze(pad:end-pad)) * gt(idx) mean(sm_aroundfreeze(pad:end-pad)) * gt(idx)], 'Color', col.stationary_sm, 'LineStyle', ':', 'LineWidth', 4)
sch = scatter(0, mean(sm_aroundfreeze(pad:end-pad)) * gt(idx), 500, 'filled', 'Clipping', 'off', 'MarkerFaceColor', col.visual_aid.mu);
set(gca ,'Layer', 'Top')

exporter(fh, paths, 'slope_sm_with_lines.pdf', 'flag', true )

%%

rng(3)
sim_params.dt = 1/pad;
sim_params.T = 30;
sim_params.time_vector = 0:sim_params.dt:sim_params.T;
sim_params.z = 0;

fh = figure('color','w', 'Position', [100 100 1100 400]);
tiledlayout(2, 3, 'TileSpacing', 'compact', 'Padding', 'compact')

% nexttile(1, [1, 2])
% hold on
% model_eval = eval(sprintf('model_%s', model));
% ncomp_vars = evaluate_model(model_eval, gt, bout_tbl);
% rt_st = nan(1,1);
% for idx_trials = 1
%     [rt_st(idx_trials), traj_st] = drift_diff_new('mu', ncomp_vars.mu2, 'theta', 0, ...
%         'z', -ncomp_vars.theta2, 'dt', sim_params.dt, 'T', sim_params.T, 'ndt',  ncomp_vars.tndt);
%     plot(traj_st, 'k-', 'LineWidth', 1)
% end
% xlim([0 30*pad])
% yticks([])
% rt_st = rt_st * pad;
% ax2 = gca;
% ax2.XAxisLocation = 'origin';
%  xticks([])
% ylabel('Safety Evidence')
% xline(0,  'k--', 'LineWidth', 2)

% apply_generic(ax2, 18);
% nexttile(1)
% histogram(rt_st, 50, 'Normalization', 'pdf')
% xlim([0 30*pad])
% ax = gca;
% yticks([])
% apply_generic(ax)

nexttile(4, [1, 2])
hold on
x_plot = (bout_tbl.onsets - pad: bout_tbl.onsets + bout_tbl.durations + pad) - bout_tbl.onsets;
plot(x_plot, sm_aroundfreeze, 'k-', 'LineWidth', 2)
plot(x_plot(pad:end-pad), sm_aroundfreeze(pad:end-pad), 'Color', 'w', 'LineWidth', 3 )
plot(x_plot(pad:end-pad), repmat(mean(sm_aroundfreeze(pad:end-pad)), length(x_plot(pad:end-pad)), 1), 'Color', col.stationary_sm, 'LineWidth', 3)

hDot_u  = plot(1, traj_st(1), 'ko', 'MarkerFaceColor','k', 'MarkerSize', 12, 'MarkerEdgeColor', col.stationary_sm, 'LineWidth', 2);  % moving dot
ax3 = gca;
ax3.XAxis.Visible = 'off';
apply_generic(ax3)
ylabel('Social Motion')
xline(0,  'k--', 'Label', 'Freeze Entry', 'LineWidth', 2, 'FontSize', 20, 'LabelOrientation', 'horizontal')
xline(bout_tbl.durations, 'k--', 'Label', 'Freeze Exit', 'LineWidth', 2, 'FontSize', 20, 'LabelOrientation', 'horizontal')
ylim([0 4])
yticks([])
apply_generic(ax3, 18);
xlim([-pad 660])
avg_tv = repmat(mean(sm_aroundfreeze(pad:end-pad)), length(x_plot(pad:end-pad)), 1);

%
nexttile(1, [1, 2])
hold on

model_eval = eval(sprintf('model_%s', model));
ncomp_vars = evaluate_model(model_eval, gt, bout_tbl);

% Run a single simulated trajectory (as you had)
[rt_sample, traj_st] = drift_diff_new('mu', ncomp_vars.mu2, 'theta', 0, ...
    'z', -ncomp_vars.theta2, 'dt', sim_params.dt, 'T', sim_params.T, 'ndt', ncomp_vars.tndt);

% Pre-create graphics objects for efficiency
hLine = plot(nan, nan, 'k-', 'LineWidth', 1);                      % growing trajectory
hDot  = plot(1, traj_st(1), 'ko', 'MarkerFaceColor','k', 'MarkerSize', 10);  % moving dot

% Choose an arrow length relative to the y-range; update after initial limits set
% (If you want a fixed length, set arrowLen to a constant.)

% Axis cosmetics (as you had)
yticks([])
ax2 = gca;
ax2.XAxisLocation = 'origin';
xticks([])
ylabel('Safety Evidence')
xlabel('Decision Bound')
xline(0,  'k--', 'LineWidth', 2)


apply_generic(ax2, 18);
linkaxes([ax2, ax3], 'x')
ylim([-5 0]);
xlim([-pad 660])

hArrow = quiver(1, traj_st(1), 0, 2*mean(sm_aroundfreeze(pad:end-pad)), ...
    'AutoScale', 'on', 'AutoScaleFactor', 1, 'MaxHeadSize', 3, ...
    'Color', col.visual_aid.mu, 'LineWidth', 1.5);

idx = find(contains(string(gt_table.Properties.VariableNames), 'mu1')); 

xxtick = [0 , 1];
nexttile(6)
hold on

plot(xxtick, xxtick * gt(idx), 'LineWidth', 3, 'Color', col.vars.sm(130,:))
ax = gca;
set(ax,'linewidth', 7,'TickDir','both');
set(ax, 'box','off')
yticks(0);
xticks([]);
xlim([0, 0.5])
ylim([-0.5, 0.5])

ax.XAxisLocation = 'origin';
apply_generic(ax, 18)
xlabel('Social Motion')
ylabel('Drift Rate (\mu)', 'Interpreter', 'tex', 'Rotation', 90)
set(ax, 'color', 'None');
set(ax,'TickLength',[0.02, 0.02])

line_v1 = plot([mean(sm_aroundfreeze(pad:end-pad)) mean(sm_aroundfreeze(pad:end-pad))], [0 mean(sm_aroundfreeze(pad:end-pad)) * gt(idx)], 'Color', col.stationary_sm, 'LineStyle', ':', 'LineWidth', 2);
line_h1 =plot([0 mean(sm_aroundfreeze(pad:end-pad))], [mean(sm_aroundfreeze(pad:end-pad)) * gt(idx) mean(sm_aroundfreeze(pad:end-pad)) * gt(idx)], 'Color', col.stationary_sm, 'LineStyle', ':', 'LineWidth', 2);
sch = scatter(0, mean(sm_aroundfreeze(pad:end-pad)) * gt(idx), 80, 'filled', 'Clipping', 'off', 'MarkerFaceColor', col.visual_aid.mu);
set(gca ,'Layer', 'Top')

% Make sure the y-limits are set before choosing arrow length
drawnow;
yl = ylim;

%
% --- Animation loop ---
n = numel(traj_st(~isnan(traj_st)));
% Optional: animate ~real-time; set to sim_params.dt for 1x speed, or smaller for faster
framePause = sim_params.dt;   % change to, e.g., 0 for fastest rendering

for k = 2:n
    % Update line to current time
    set(hLine, 'XData', 1:k, 'YData', traj_st(1:k));

    % Move dot
    set(hDot, 'XData', k, 'YData', traj_st(k));

    set(hDot_u, 'XData', k, 'YData', avg_tv(k));

    % Keep arrow at the dot, pointing up
    set(hArrow, 'XData', k, 'YData', traj_st(k), 'VData', mean(sm_aroundfreeze(pad:end-pad)));

    % Render
    drawnow limitrate
    if framePause > 0
        pause(framePause);
    end
end

% Convert RT to seconds if you still need it later
rt_st = rt_sample * pad;


