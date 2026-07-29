clearvars

% Mean social motion around the start of a loom-evoked freezing bout.
% Offset-aligned counterpart: sm_aligned_to_freeze_offset.
%
% t = 0 is the bout onset. Each panel carries two traces, and the gap between
% them is the point:
%
%   sm_freeze  only the frames the fly was actually immobile, so at lag k the
%              mean is taken over the bouts still running at k. The sample
%              shrinks as the lag grows
%   sm_ili     a fixed window from the onset, ignoring where the bout ended.
%              Every bout contributes at every lag, so the sample is constant
%              and the later part mixes frozen and post-freeze frames
%
% Reading them together separates a real change in social motion during a
% freeze from one that only reflects which bouts are still running.
%
% CENSORING. The onset anchor is always genuine, so unlike the offset script
% censoring does not corrupt t = 0. What it does is cap how far forward each
% bout's data can be trusted, which is what ending_time records. So every bout
% is kept and its during-freeze trace is truncated at its own ending_time,
% rather than whole bouts being dropped.
%
% Window-censored bouts are deliberately kept here even though the offset
% script drops them. Dropping them would remove exactly the longest freezes,
% which is the survivorship artefact this script exists to expose; and their
% frames up to the window edge are perfectly good observations. They are split
% into their own panel so the truncation stays visible rather than being
% averaged in silently.

% Load the table first. We will take advantage of an already existing
% dataset.
threshold_imm = 2; threshold_mob = 2; threshold_pc = 4; id_code = sprintf('imm%d_mob%d_pc%d', threshold_imm, threshold_mob, threshold_pc);
paths = path_generator('folder', 'social_motion/proximity', 'bouts_id', id_code, 'imfirst', false);
bouts = importdata(fullfile(paths.dataset, 'bouts.mat'));

% Here you can change the window of loom-evoked
thresholds = define_thresholds('le_window', struct('le_window_sl', [5 55], 'le_window_fl', [5 55]));
%thresholds = define_thresholds;

bouts = bouts_formatting(bouts, thresholds);

% Process the data: set thresholds for sm, fs and ln. Set minimum duration.
bouts_proc = data_parser_new(bouts, 'type', 'immobility', 'period', 'loom', 'window', 'le', 'nloom', 2:20, 'min_dur', 30);

% Adds is_censored / censored_contacts / censored_loom / ending_time. Slow: it
% loads the mindist cache and walks every bout, so call it once
bouts_proc = impose_contact_threshold(bouts_proc);

inizio = 0;
fine = 630;

uncensored = ~bouts_proc.is_censored;
collision  = bouts_proc.censored_contacts;
window_cut = bouts_proc.censored_loom;

fprintf('%d bouts: %d uncensored, %d collision-censored, %d window-censored\n', ...
    height(bouts_proc), sum(uncensored), sum(collision), sum(window_cut));

% Retargeting ends to onsets + ending_time truncates each bout at the point it
% stopped being observed cleanly: the contact frame for a collision, the window
% edge for a window-censored bout. ending_time equals durations for uncensored
% bouts, so this leaves those untouched
bouts_trunc = bouts_proc;
bouts_trunc.ends = bouts_proc.onsets + bouts_proc.ending_time;

panel_mask = {uncensored, collision, window_cut};
panel_name = {'Uncensored', 'Collision-censored', 'Window-censored'};

t_ili = inizio : fine;

fh = figure('color', 'w', 'Position', [100 100 1500 450]);
tlo = tiledlayout(1, numel(panel_mask), 'TileSpacing', 'compact', 'Padding', 'compact');
ax = gobjects(1, numel(panel_mask));

for p = 1:numel(panel_mask)

    sel = panel_mask{p};
    ax(p) = nexttile;
    hold on

    % extract_sm_from_bouts builds nan(0, []) on an empty table and errors, so
    % an empty selection has to be caught before it gets there
    if ~any(sel)
        title(sprintf('%s (n = 0)', panel_name{p}));
        continue
    end

    b = bouts_trunc(sel, :);

    % Truncated at ending_time by the retargeting above, so this trace never
    % runs past the point the bout was still cleanly observed
    sm_freeze = extract_sm_from_bouts(b, 'type', 'onlyfreeze', 'output_type', 'mat', 'align', 'onset');
    sm_freeze = sm_freeze(:, 1:min(end, fine));

    % Left untruncated on purpose: this is the boundary-ignoring control, and
    % capping it at ending_time would make it a copy of the trace above
    sm_ili = extract_sm_from_bouts(b, 'type', 'onsets', 'output_type', 'mat', 'window', [inizio fine]);

    t_freeze = 0 : size(sm_freeze, 2) - 1;

    plot(t_freeze, mean(sm_freeze, 1, 'omitnan'), 'k-', 'LineWidth', 2)
    plot(t_ili, mean(sm_ili, 1, 'omitnan'), 'r-', 'LineWidth', 2)

    xline(0, '--k', 'HandleVisibility', 'off');
    xlabel('Frames from freeze onset')
    ylabel('Social motion')
    title(sprintf('%s (n = %d)', panel_name{p}, sum(sel)));

    if p == 1
        legend({'during freeze only', 'fixed window'}, 'Location', 'northeast')
    end
end

linkaxes(ax, 'xy')
