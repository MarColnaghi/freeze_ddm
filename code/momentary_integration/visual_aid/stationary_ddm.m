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

    lh.Color = col.timevarying_sm;
    figure(fh)
    plot(x_plot(pad:end-pad), repmat(mean(sm_aroundfreeze(pad:end-pad)), length(x_plot(pad:end-pad)), 1), 'Color', col.stationary_sm, 'LineWidth', 3)
    %yticks(sort([ax.YTick round(mean(sm_aroundfreeze(pad:end-pad)),2)]))
    exporter(fh, paths, 'st_and_tv_signal.pdf', 'flag', true )

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
sch2 = scatter(0, mean(sm_aroundfreeze(pad:end-pad)) * gt(idx), 500, 'filled', 'Clipping', 'off', 'MarkerFaceColor', col.visual_aid.mu);
set(gca ,'Layer', 'Top')

exporter(fh, paths, 'slope_sm_with_lines.pdf', 'flag', true )




%% --- Setup ----------------------------------------------------------------

cd(paths.fig)
rng(56);
tv = false;


assert(exist('pad','var')==1 && pad>0, 'pad must be a positive scalar.');
assert(exist('model','var')==1, 'model must exist (string).');
assert(exist('gt','var')==1, 'gt must exist.');
assert(exist('bout_tbl','var')==1 && istable(bout_tbl), 'bout_tbl must be a table.');
assert(exist('sm_aroundfreeze','var')==1, 'sm_aroundfreeze must exist.');
assert(exist('evaluate_model','file')==2, 'evaluate_model.m not found.');
assert(exist('drift_diff_new','file')==2, 'drift_diff_new.m not found.');
assert(exist('apply_generic','file')==2, 'apply_generic.m not found.');

sim_params.dt = 1/60;
sim_params.T  = 30;
sim_params.time_vector = 0:sim_params.dt:sim_params.T;
sim_params.z  = 0;

% use first gt entry as slope/beta for mapping panel
beta = gt(1);

% model function handle (avoids eval)
model_fn   = str2func(['model_' char(model)]);
mo = model_fn();
ncomp_vars = evaluate_model(mo, gt, bout_tbl);

% --- Simulate a single trajectory ----------------------------------------

tv_signal = motion_cache(bout_tbl.fly);
tv_signal = tv_signal(bout_tbl.onsets:bout_tbl.onsets + length(sim_params.time_vector) - 2) .* beta;


for idx_tries = 1:1000
    
    [rt_sample_tv(idx_tries)] = drift_diff_new('mu_t', tv_signal, 'theta', 0, ...
        'z', -ncomp_vars.theta2, 'dt', sim_params.dt, 'T', sim_params.T, 'ndt', 0, 'seed', idx_tries);
    [rt_sample_st(idx_tries)] = drift_diff_new('mu', ncomp_vars.mu2, 'theta', 0, ...
        'z', -ncomp_vars.theta2, 'dt', sim_params.dt, 'T', sim_params.T, 'ndt', 0, 'seed', idx_tries);
end

[d_tv, i_tv] = min(abs(rt_sample_tv - bout_tbl.durations ./60));
[d_st, i_st] = min(abs(rt_sample_st - bout_tbl.durations ./60));

if tv
    [rt_sample_tv, traj] = drift_diff_new('mu_t', tv_signal, 'theta', 0, ...
        'z', -ncomp_vars.theta2, 'dt', sim_params.dt, 'T', sim_params.T, 'ndt', ncomp_vars.tndt, 'seed', i_tv);
else
    [rt_sample_st, traj] = drift_diff_new('mu', ncomp_vars.mu2, 'theta', 0, ...
        'z', -ncomp_vars.theta2, 'dt', sim_params.dt, 'T', sim_params.T, 'ndt', ncomp_vars.tndt, 'seed', i_st);
end

traj_valid = traj(~isnan(traj));
n_traj     = numel(traj_valid);

% --- Build x over freeze window ------------------------------------------
x_plot = (bout_tbl.onsets - pad : bout_tbl.onsets + bout_tbl.durations + pad) - bout_tbl.onsets;
% guard in case pad larger than trace ends
core_idx = (pad+1) : (numel(x_plot)-pad);
core_idx = core_idx(core_idx>=1 & core_idx<=numel(sm_aroundfreeze)); % clamp
avg_core = mean(sm_aroundfreeze(core_idx));

% stationary vs time-varying SM signals (match animation length later)
sm_signal_tv = sm_aroundfreeze(core_idx);
sm_signal_st = ones(size(sm_signal_tv)) * avg_core;

try
    vw = VideoWriter('st.mp4','MPEG-4');   % .mp4 (H.264)
catch
    vw = VideoWriter('st.avi','Motion JPEG AVI');  % fallback
    vw.Quality = 95;  % Only for Motion JPEG AVI
end

vw.FrameRate = round(1/sim_params.dt);  % ~60 fps from your dt=1/60
open(vw);

% --- Figure & layout ------------------------------------------------------
fh = figure('color','w', 'Position', [100 100 1300 500]);
tiledlayout(2, 3, 'TileSpacing', 'loose', 'Padding', 'loose');

% --- (A) Social Motion panel ---------------------------------------------
ax_sm = nexttile(2, [1, 2]); hold(ax_sm,'on');

plot(ax_sm, x_plot, sm_aroundfreeze, 'k-', 'LineWidth', 2);
plot(ax_sm, x_plot(core_idx), sm_aroundfreeze(core_idx), 'Color','w', 'LineWidth', 3);

if tv
    plot(ax_sm, x_plot(core_idx), sm_signal_tv, 'Color', col.timevarying_sm, 'LineWidth', 3);
    hDot_sm = plot(ax_sm, x_plot(core_idx(1)), sm_aroundfreeze(core_idx(1)), 'ko', ...
        'MarkerFaceColor', 'w', 'MarkerSize', 12, 'MarkerEdgeColor', col.timevarying_sm, 'LineWidth', 2);
else
    plot(ax_sm, x_plot(core_idx), sm_signal_st, 'Color', col.stationary_sm, 'LineWidth', 3);
    hDot_sm = plot(ax_sm, x_plot(core_idx(1)), sm_aroundfreeze(core_idx(1)), 'ko', ...
        'MarkerFaceColor', 'w', 'MarkerSize', 12, 'MarkerEdgeColor', col.stationary_sm, 'LineWidth', 2);
end

ax_sm.XAxis.Visible = 'off';
apply_generic(ax_sm);
ylabel(ax_sm,'Social Motion');
xline(ax_sm, 0, 'k--', 'Label','Freeze Entry', 'LineWidth',2, 'FontSize',20, 'LabelOrientation','horizontal');
xline(ax_sm, bout_tbl.durations(1), 'k--', 'Label','Freeze Exit', 'LineWidth',2, 'FontSize',20, 'LabelOrientation','horizontal');
yticks(ax_sm, []);
apply_generic(ax_sm, 24);
set(ax_sm,'linewidth', 4, 'TickDir','out');

ylim(ax_sm, [0 4]); xlim(ax_sm, [-pad 660]);

% --- (B) Drift mapping panel ---------------------------------------------
ax_map = nexttile(1, [1 1]); hold(ax_map,'on');

plot(ax_map, [0 2], [0 2]*beta, 'LineWidth', 4, 'Color', col.vars.sm(130,:));

set(ax_map,'LineWidth',1.5,'TickDir','both','Box','off','Color','none','TickLength',[0.02 0.02]);
yticks(ax_map, 0); xticks(ax_map, []);
xlim(ax_map, [0 2]); ylim(ax_map, [-1.5 1.5]);
ax_map.XAxisLocation = 'origin';
apply_generic(ax_map, 24);
set(ax_map,'linewidth', 4, 'TickDir','out');
xlabel(ax_map,'Social Motion');
ylabel(ax_map,'Drift Rate (\mu)', 'Interpreter','tex','Rotation',90);

% Guides at mean SM
mx = avg_core; my = beta*avg_core;
line_v1 = plot(ax_map, [mx mx], [0 my], 'Color', col.stationary_sm, 'LineStyle',':', 'LineWidth', 2);
line_h1 = plot(ax_map, [0 mx], [my my], 'Color', col.stationary_sm, 'LineStyle',':', 'LineWidth', 2);
if tv
    sch1     = scatter(ax_map, mx, - 0.1, 250, 'filled', 'Clipping','off', 'MarkerFaceColor', col.timevarying_sm, 'Marker', '^');
else
    sch1     = scatter(ax_map, mx, - 0.1, 250, 'filled', 'Clipping','off', 'MarkerFaceColor', col.stationary_sm, 'Marker', '^');
end
sch2     = scatter(ax_map, 0, my, 200, 'filled', 'Clipping','off', 'MarkerFaceColor', col.visual_aid.mu);
set(ax_map,'Layer','Top');

% --- (C) DDM panel --------------------------------------------------------
ax_ddm = nexttile(5, [1, 2]); hold(ax_ddm,'on');

hLine = plot(ax_ddm, nan, nan, 'k-', 'LineWidth', 1.5);
hDot  = plot(ax_ddm, 1, traj(1), 'ko', 'MarkerFaceColor','k', 'MarkerSize', 14);

yticks(ax_ddm, []);
ax_ddm.XAxisLocation = 'origin';
xticks(ax_ddm, []);
ylabel(ax_ddm,'Safety Evidence');
xlabel(ax_ddm,'Decision Bound');
xline(ax_ddm, 0, 'k--', 'LineWidth', 2);
apply_generic(ax_ddm, 24);
set(ax_ddm,'linewidth', 4, 'TickDir','out');
ylim(ax_ddm, [-5 0]); xlim(ax_ddm, [-pad 660]);

hArrow = quiver(ax_ddm, 1, traj(1), 0, avg_core, ...
    'AutoScale','on','AutoScaleFactor',1,'MaxHeadSize',3, ...
    'Color', col.visual_aid.mu, 'LineWidth', 2, 'Clipping', 'off');

linkaxes([ax_ddm, ax_sm], 'x');
drawnow; yl = ylim(ax_ddm);

% --- Animation ------------------------------------------------------------
n_sm = numel(sm_signal_st);
n    = n_traj;
framePause = sim_params.dt;  % 0 for fastest

if tv
    signal = sm_signal_tv;
else
    signal = sm_signal_st;
end

for k = 2:n
    % DDM
    set(hLine, 'XData', 1:k, 'YData', traj(1:k));
    set(hDot,  'XData', k,   'YData', traj(k));

    % Clamp SM index to available samples
    kk = min(k, n_sm);
    set(hDot_sm, 'XData', x_plot(core_idx(kk)), 'YData', signal(kk));

    % Update guide point/lines in mapping panel using current stationary SM
    cur_sm = signal(kk);
    cur_mu = beta * cur_sm;
    set(line_v1, 'XData', [cur_sm cur_sm], 'YData', [0 cur_mu]);
    set(line_h1, 'XData', [0 cur_sm],      'YData', [cur_mu cur_mu]);
    set(sch1,     'XData', cur_sm,          'YData', -0.1);
    set(sch2,     'XData', 0,               'YData', cur_mu);

    % Arrow at DDM dot with stationary SM magnitude
    set(hArrow, 'XData', k, 'YData', traj(k), 'VData', cur_mu * 3);

    drawnow limitrate;
    if framePause>0, pause(framePause); end

    if k > 270 && k < 330
        pause(0)

    elseif k == n
        set(hArrow, 'Visible', 'off');

    end

    writeVideo(vw, getframe(fh));

end
close(vw);


% Output RT in seconds (if pad is frames/sec)
%rt_st_seconds = rt_sample * pad; 
%%

fh = figure('color','w', 'Position', [100 100 900 500]);
tiledlayout(3, 1, 'TileSpacing', 'none', 'Padding', 'loose');

nth = nexttile(1);

nexttile(2)
ax_ddm = gca;

hold on
yticks(ax_ddm, []);
ax_ddm.XAxisLocation = 'origin';
xticks(ax_ddm, []);
ylabel(ax_ddm,'Safety Evidence');
xline(ax_ddm, 0, 'k--', 'LineWidth', 2);
apply_generic(ax_ddm, 18);
set(ax_ddm,'linewidth', 4, 'TickDir','out');
ylim(ax_ddm, [-ncomp_vars.theta2*2 0]); xlim(ax_ddm, [-1 15]);

for idx_tries = 1:20000
    
    [rt_sample_tv(idx_tries), traj_tv] = drift_diff_new('mu_t', tv_signal, 'theta', 0, ...
        'z', -ncomp_vars.theta2, 'dt', sim_params.dt, 'T', sim_params.T, 'ndt', 0, 'seed', idx_tries);
    [rt_sample_st(idx_tries), traj_st] = drift_diff_new('mu', ncomp_vars.mu2, 'theta', 0, ...
        'z', -ncomp_vars.theta2, 'dt', sim_params.dt, 'T', sim_params.T, 'ndt', 0, 'seed', idx_tries);

    if idx_tries < 500
        plot(sim_params.time_vector,traj_tv, 'Color', [hex2rgb(col.timevarying_sm) 0.05])
        plot(sim_params.time_vector, traj_st, 'Color', [hex2rgb(col.stationary_sm) 0.05])
    end

end
scatter(0, -ncomp_vars.theta2, 80, 'k', 'filled')

nexttile(1)
hold on
histogram(rt_sample_tv, 1/120:1/60*6:15, 'FaceColor', col.timevarying_sm, 'Normalization', 'pdf', 'EdgeColor', 'none')
histogram(rt_sample_st,  1/120:1/60*6:15, 'FaceColor', col.stationary_sm, 'Normalization', 'pdf', 'EdgeColor', 'none')
ax = gca;
ax.XAxis.Visible = 'on'; xticks([]);
ax.YAxis.Visible = 'off';

set(ax ,'Layer', 'Top')

nexttile(3)
hold on
d = motion_cache(bout_tbl.fly);
c = d(bout_tbl.onsets:bout_tbl.onsets + length(sim_params.time_vector) - 2);
x_plot = sim_params.time_vector(1:length(c ));
lh = plot(x_plot, c, 'Color', 'k', 'LineWidth', 2);
ax3 = gca;
ax3.XAxis.Visible = 'off';
ax3.YAxis.Visible = 'on'; ax3.YAxisLocation = 'right';
apply_generic(ax3, 18)
ylabel('Social Motion'); yticks([0 5]);
ylim([-2 5])
linkaxes([ax_ddm, ax, ax3], 'x')
xlim(ax_ddm, [-1 15]);

exporter(fh, paths, 'simulation_example.pdf', 'flag', true )
