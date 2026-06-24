%% freeze_break_sm_grand_average.m
% -------------------------------------------------------------------------
% Grand-average of the social motion the focal fly experiences, aligned to
% the moment it BREAKS its freeze (un-freeze / freeze offset).
%
% One clean panel meant for a presentation slide: bold grand-average line
% (mean across freezes) with a shaded SEM band, time measured in seconds
% relative to the freeze break (t = 0).
%
% Built from the same primitives as signals_aligned_2_offset.m, but trimmed
% down to a single, slide-friendly figure.
%
% Notes / analysis decisions (see project memory):
%   - exclude loom 1, use nloom 2:20
%   - loom-evoked freezes only, min duration 0.5 s (30 frames)
%   - right-align each freeze to its break, then average column-wise
%   - drop collision-contaminated freezes (min distance < contact threshold)
%     because contact spikes pixel-change and corrupts the signal near the
%     break, which is exactly the epoch we care about here.
% -------------------------------------------------------------------------

clearvars

%% ---- Parameters ---------------------------------------------------------
camera_fps      = 60;
delay_after     = 90;     % frames of signal to keep AFTER the break (1.5 s)
contact_thresh  = 80;     % px; conservative collision threshold (censoring)
max_dur         = 630;    % frames; drop freezes longer than the obs. window
min_n           = 30;     % only plot lags with at least this many freezes
xlims_s         = [-6 delay_after/camera_fps];   % seconds, relative to break
sm_cache        = 'motion_cache';                % social-motion timeseries

%% ---- Load & parse bouts -------------------------------------------------
threshold_imm = 2; threshold_mob = 2; threshold_pc = 4;
id_code = sprintf('imm%d_mob%d_pc%d', threshold_imm, threshold_mob, threshold_pc);
paths = path_generator('folder', 'social_motion', 'bouts_id', id_code, 'imfirst', false);

bouts = importdata(fullfile(paths.dataset, 'bouts.mat'));

% Loom-evoked window (mirrors signals_aligned_2_offset.m)
thresholds = define_thresholds('le_window', ...
    struct('le_window_sl', [5 55], 'le_window_fl', [5 55]));
bouts = bouts_formatting(bouts, thresholds);

% Loom-evoked freezes, exclude loom 1, min duration 0.5 s
bouts_proc = data_parser_new(bouts, 'type', 'immobility', 'period', 'loom', ...
    'window', 'le', 'nloom', 2:20, 'min_dur', 30);

%% ---- Extract break-aligned social motion --------------------------------
% Right-aligned matrix: each row is one freeze, last columns sit at the break
% (+ delay_after). NaNs pad the front of shorter freezes.
sm_freeze = extract_sm_from_bouts(bouts_proc, 'type', 'onlyfreeze', ...
    'output_type', 'mat', 'align', 'offset', 'delay_after', delay_after, ...
    'cache', sm_cache);

%% ---- Censor collision-contaminated freezes ------------------------------
[~, contact_mask] = impose_contact_threshold(bouts_proc, 'threshold', contact_thresh);
keep_mask = ~contact_mask & (bouts_proc.durations <= max_dur);
sm_keep   = sm_freeze(keep_mask, :);

%% ---- Grand average ------------------------------------------------------
% Time axis in frames relative to the break (t = 0), then seconds.
t_frames = (-size(sm_keep, 2) + 1 : 1 : 0) + delay_after;
t_s      = t_frames / camera_fps;

n_contrib = sum(~isnan(sm_keep), 1);
avg = mean(sm_keep, 1, 'omitnan');
sem = std(sm_keep, 0, 1, 'omitnan') ./ sqrt(n_contrib);

% Only show lags supported by enough freezes (avoids a noisy deep-past tail)
valid = n_contrib >= min_n & t_s >= xlims_s(1) & t_s <= xlims_s(2);
t_s   = t_s(valid);
avg   = avg(valid);
sem   = sem(valid);

%% ---- Plot ---------------------------------------------------------------
col  = cmapper;
line_col = col.var.sm;          % social-motion red

fh = figure('color', 'w', 'Position', [100 100 760 520]);
ax = axes(fh); hold(ax, 'on');

% SEM band
fill([t_s, fliplr(t_s)], [avg + sem, fliplr(avg - sem)], line_col, ...
    'EdgeColor', 'none', 'FaceAlpha', 0.20, 'HandleVisibility', 'off');

% Grand-average line
plot(t_s, avg, 'LineWidth', 3, 'Color', line_col);

% Break marker
yl = ylim(ax);
xline(0, '--k', 'LineWidth', 1.5, 'HandleVisibility', 'off');

apply_generic(ax, 'font_size', 20, 'tick_length', 0.015, 'line_width', 2, ...
    'xlims', xlims_s, 'xticks', xlims_s(1):1:xlims_s(2));

xlabel('Time from freeze break (s)')
ylabel('Social motion')
title(sprintf('Break-aligned social motion (n = %d freezes)', sum(keep_mask)), ...
    'FontWeight', 'normal', 'FontSize', 18)

% Label the break tick
ticks  = ax.XAxis.TickValues;
labels = arrayfun(@(v) num2str(v), ticks, 'UniformOutput', false);
labels(ticks == 0) = {'Break'};
xticklabels(ax, labels)

%% ---- Export -------------------------------------------------------------
exporter(fh, paths, 'freeze_break_sm_grand_average.pdf')
% PNG for slides (raster — exportgraphics directly, exporter forces vector)
exportgraphics(fh, fullfile(paths.fig, 'freeze_break_sm_grand_average.png'), ...
    'BackgroundColor', 'white', 'Resolution', 300);
