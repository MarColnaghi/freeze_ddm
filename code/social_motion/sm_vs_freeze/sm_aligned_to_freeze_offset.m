clearvars

% Mean social motion around the end of a loom-evoked freezing bout.
% Onset-aligned counterpart: sm_aligned_to_freeze_onset.
%
% t = 0 is the last observed frozen frame, so negative lags are inside the
% freeze and positive lags are after it. Each panel carries two traces, and the
% gap between them is the point:
%
%   sm_freeze   only the frames the fly was actually immobile, right-aligned on
%               t = 0. At lag -k the mean is taken over the bouts that lasted at
%               least k frames, so the sample shrinks the further back you look
%   sm_window   a fixed window ending at t = 0, ignoring where the bout began.
%               Every bout contributes at every lag, so the sample is constant
%               and the early part mixes pre-freeze and frozen frames
%
% Reading them together separates a real change in social motion before a
% freeze breaks from one that only reflects which bouts survived that long.
%
% CENSORING. Offset alignment is far more sensitive to it than onset alignment,
% because censoring lands on the anchor itself rather than on the tail. The
% three populations are not three views of one thing:
%
%   uncensored        t = 0 is a genuine, spontaneous freeze break
%   collision         t = 0 is the frame another fly came within threshold. A
%                     real event, but a different question
%   window-censored   the bout outlasted the observation window, so no event
%                     exists to align on. Counted and dropped, since aligning
%                     them would put an arbitrary timestamp at t = 0
%
% Pooling the last two, as the duration histograms reasonably do, would put a
% real event and a non-event in the same average.

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

fine = 630;
delay_after = 90;   % frames kept past the anchor, as in signals_aligned_2_offset

uncensored = ~bouts_proc.is_censored;
collision  = bouts_proc.censored_contacts;
window_cut = bouts_proc.censored_loom;

fprintf('%d bouts: %d uncensored, %d collision-censored, %d window-censored (dropped)\n', ...
    height(bouts_proc), sum(uncensored), sum(collision), sum(window_cut));

% For the collision panel the anchor is the contact frame, not the recorded end
% of the bout. ending_time is the 1-based frame within the bout at which contact
% happened, so retargeting ends to onsets + ending_time makes 'align', 'offset'
% put t = 0 there. For uncensored bouts ending_time equals durations, so the
% same expression leaves that panel's anchor exactly where it was
bouts_at_contact = bouts_proc;
bouts_at_contact.ends = bouts_proc.onsets + bouts_proc.ending_time;

panel_mask  = {uncensored, collision};
panel_table = {bouts_proc, bouts_at_contact};
panel_name  = {'Uncensored (freeze break)', 'Collision (contact)'};

n_keep = fine + delay_after + 1;
t_window = -fine : delay_after;

fh = figure('color', 'w', 'Position', [100 100 1100 450]);
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

    b = panel_table{p}(sel, :);

    sm_freeze = extract_sm_from_bouts(b, 'type', 'onlyfreeze', 'output_type', 'mat', ...
        'align', 'offset', 'delay_after', delay_after);

    % Trim to the common window. Anchoring on the right edge keeps t = 0 on the
    % anchor whatever the longest bout happens to be, and the max(1, ...) leaves
    % the matrix alone when no bout is long enough to fill the window
    sm_freeze = sm_freeze(:, max(1, end - n_keep + 1):end);

    % extract_sm_from_bouts windows the 'onsets' type off bouts.onsets and has
    % no offset equivalent, but that branch reads only .onsets and .fly and
    % NaN-fills any index outside the recording, so retargeting the column is
    % enough to move the window onto the anchor. ends is one past the anchor
    b_win = b;
    b_win.onsets = b.ends - 1;
    sm_window = extract_sm_from_bouts(b_win, 'type', 'onsets', 'output_type', 'mat', ...
        'window', [-fine delay_after]);

    t_freeze = (-size(sm_freeze, 2) + 1 : 0) + delay_after;

    plot(t_freeze, mean(sm_freeze, 1, 'omitnan'), 'k-', 'LineWidth', 2)
    plot(t_window, mean(sm_window, 1, 'omitnan'), 'r-', 'LineWidth', 2)

    xline(0, '--k', 'HandleVisibility', 'off');
    xlabel('Frames from anchor')
    ylabel('Social motion')
    title(sprintf('%s (n = %d)', panel_name{p}, sum(sel)));

    if p == 1
        legend({'during freeze only', 'fixed window'}, 'Location', 'northwest')
    end
end

linkaxes(ax, 'xy')
