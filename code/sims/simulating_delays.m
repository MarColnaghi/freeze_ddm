%% ════════════════════════════════════════════════════════════════════════
%  SIMULATING DELAYS — what each delay parameter actually does
%  ------------------------------------------------------------------------
%  The model has delays that are easy to confuse because they all push reaction
%  times later. They act through completely different mechanisms:
%
%    sensory_delay  -> shifts the INPUT SIGNAL (evidence arrives late).
%                      Changes the drive.
%    delayed_start  -> leaves the signal untouched but delays WHEN THE
%                      ACCUMULATOR STARTS integrating. Changes the path.
%    ndt            -> touches neither signal nor path; a pure POST-HOC shift
%                      of the reported RT. Changes only the readout.
%
%  ...and one hidden ASSUMPTION, which gets its own column:
%
%    With sensory_delay = L, the evidence at freeze-relative time tau is
%    raw[tau - L]. For tau < L that is BEFORE freeze onset. bayes_fpe fills that
%    window with ZEROS -- not physiology, just an artefact of its signal array
%    starting at onset. So zero-padding makes sensory_delay do DOUBLE DUTY: it
%    delays evidence AND deletes the first L seconds of it. The fly's visual
%    stream did not begin at freeze onset, and that motion is in the cache, so we
%    can feed the REAL pre-onset stream instead (p.n_pre) and make sensory_delay
%    a PURE lag. Columns 1 and 2 differ ONLY in that choice.
%
%  Layout: one column per condition, three rows showing where each becomes visible.
%    Row 1 (input drift)     -> only the sensory columns move it (and they differ
%                               from each other on tau in [0, L)).
%    Row 2 (trajectories)    -> sensory + delayed_start move it.
%    Row 3 (RT distribution) -> all of them move it.
%  Trajectories use the SAME SEEDS across every condition, so any difference is
%  the manipulation, not the noise.
%
%  Driven by a real freeze's social-motion trace, loaded with a pre-onset margin.
% ════════════════════════════════════════════════════════════════════════
clearvars; clc;

this_dir = fileparts(mfilename('fullpath'));
addpath(fullfile(this_dir, 'simulators'));
addpath(fullfile(this_dir, 'simulators', 'cpp_mex_code'));

% ── Model regime ─────────────────────────────────────────────────────────
dt     = 1/60;          % s / frame (native social-motion rate)
sigma  = 1;           % diffusion SD
theta  = 2.7;           % bound -> median RT ~1.1 s on this trial (delays visible)
beta   = 4.2;           % drift gain (rate per unit social motion)
lambda = 0;             % leak (0 = perfect integrator)
ivf    = 0;           % start point as a FRACTION of theta (bayes_fpe convention)
x0     = 0;   % ...and the same start point in ABSOLUTE units.
                        % simulate_freeze_ddm takes the fraction (and converts
                        % internally); sim_leaky_accumulator/sim_ddm_seeded take
                        % the absolute value. Keep both in sync from this one knob.

idx_trial   = 5422;412;6000;%5350;%1223;     % a real freeze with a healthy motion trace
N           = 2000000;   % trials for the RT distributions
seeds_traj  = [74 211 300];   % example trajectories (SAME seeds in every condition)

delay_vals = [0 0.4 0.8];                              % s
cols = [0.20 0.20 0.20; 0.90 0.50 0.13; 0.20 0.45 0.75];

% One column per condition. The two sensory columns run IDENTICAL parameters and
% differ only in use_pre, which isolates the padding assumption as the single
% manipulated variable.
conds = struct( ...
    'type',    {'sensory_delay',            'sensory_delay',                 'delayed_start',                 'ndt'}, ...
    'use_pre', {false,                      true,                            false,                           false}, ...
    'name',    {'sensory delay (zero-pad)', 'sensory delay (real pre-onset)', 'delayed start of accumulation', 'non-decision time'});
nC = numel(conds);

% ── Real social-motion drive (one freeze), WITH a pre-onset margin ───────
% Margin matches the largest delay (+ a guard for the fractional-frame interp).
n_pre     = round(max(delay_vals)/dt) + 2;
win_start = -n_pre;                      % signed FRAME offset relative to onset
win_end   = 630;

paths      = path_generator('folder', 'sims/sanity_checks');
bouts      = importdata(fullfile(paths.dataset, 'bouts.mat'));
bouts_proc = data_parser_new(bouts, 'type','immobility', 'period','loom', ...
                             'window','le', 'nloom', 2:20);
bouts_proc.smp       = ones(height(bouts_proc), 1);
bouts_proc.intercept = ones(height(bouts_proc), 1);
% window is [start end] in FRAMES relative to the onset; a negative start returns
% pre-onset samples (same idiom as data_parser_new.m:116). Column j holds frame
% offset win_start+(j-1), so the onset sample sits at column 1-win_start.
sm_signal  = extract_sm_from_bouts(bouts_proc, 'type','onsets', ...
                             'output_type','mat', 'window', [win_start win_end], ...
                             'norm_factor', 10);

sm_raw = sm_signal(idx_trial, :).';
i0     = 1 - win_start;                  % column of the onset sample (t = 0)

% Out-of-range samples come back as NaN (extract_sm_from_bouts.m:79-87). Note
% fillmissing('previous') CANNOT fill a LEADING NaN -- left alone it would end up
% as 0 and be indistinguishable from genuine zero social motion, which is exactly
% the artefact this figure is about. So fill both directions and then assert.
n_nan_pre = sum(isnan(sm_raw(1:i0-1)));
sm_full   = fillmissing(sm_raw, 'previous');
sm_full   = fillmissing(sm_full, 'next');
assert(~any(isnan(sm_full)), 'simulating_delays: trial %d has unfillable gaps.', idx_trial);

sm_post = sm_full(i0:end);               % post-onset only (== the old [0 630] call)
n_fr    = numel(sm_post);
t_sig   = (0:n_fr-1)' * dt;              % freeze-relative time, t=0 at freeze onset
t_full  = ((0:numel(sm_full)-1)' - n_pre) * dt;   % signed axis, negative = pre-onset
assert(abs(t_full(i0)) < 1e-12, 'simulating_delays: onset column mis-aligned.');

base = struct('dt', dt, 'theta', theta, 'sigma', sigma, 'leak', lambda, ...
              'gain', beta, 'bias', 0, 'backend', 'mex', 'seed', 42, 'x0', x0);

xmax = 10.5;   % s, plotting window

fprintf('Trial %d | theta=%.2f beta=%.2f | %d post-onset frames (%.1f s)\n', ...
        idx_trial, theta, beta, n_fr, n_fr*dt);
fprintf('Pre-onset margin: %d frames (%.2f s) | NaNs in pre-onset window: %d\n', ...
        n_pre, n_pre*dt, n_nan_pre);
% The realistic contaminant is not the recording boundary but the pre-onset window
% reaching back into the PREVIOUS freeze -- report it rather than assume.
if ismember('time_since_last', bouts_proc.Properties.VariableNames)
    tsl = abs(double(bouts_proc.time_since_last(idx_trial)));
    fprintf('time_since_last = %g frames (margin %d) -> pre-onset window %s\n', ...
        tsl, n_pre, ternary(tsl >= n_pre, 'clear of the previous bout', ...
                            'MAY OVERLAP the previous bout'));
end

%% ── Figure: 4 conditions x 3 views ──────────────────────────────────────
figure('Color','w','Position',[30 30 1900 950]);
tiledlayout(3, nC, 'TileSpacing','compact', 'Padding','compact');

RT = cell(nC, numel(delay_vals));

sig_ylim  = [0 12];
traj_ylim = [-2 theta+0.5];
pdf_ylim  = [0 2];
% subtitle() is placed by the layout itself, so the note can never drift onto a
% neighbouring panel the way a manually positioned text() can.
note = @(ax, str) subtitle(ax, str, 'FontSize', 12, 'FontAngle','italic', ...
                           'Color', [0.35 0.35 0.35]);

% Each column builds its params inline (2 lines, repeated per row): the zero-pad
% columns get the post-onset signal only (n_pre = 0, i.e. the bayes_fpe view of the
% world); the pre-onset column gets the full trace plus n_pre so the delay can
% reach back into real evidence.
%
% Row 1 spans the pre-onset margin so the evidence source is visible.
sig_xlim = [-n_pre*dt xmax];

% ---------- Row 1: the input drift the accumulator actually sees ----------
for c = 1:nC
    ax = nexttile(c); hold(ax, 'on')
    % shade t < 0: the pre-onset interval, where the two sensory columns disagree
    patch(ax, [sig_xlim(1) 0 0 sig_xlim(1)], [sig_ylim(1) sig_ylim(1) sig_ylim(2) sig_ylim(2)], ...
          [0.93 0.93 0.93], 'EdgeColor','none', 'HandleVisibility','off');
    for v = 1:numel(delay_vals)
        p = base; p.(conds(c).type) = delay_vals(v);
        if conds(c).use_pre, sig = sm_full; p.n_pre = n_pre;
        else,                sig = sm_post; p.n_pre = 0; end
        [drift, full_d] = signal_to_drift(sig, p);
        % plot the untrimmed drift so the pre-onset region shows where the
        % evidence actually comes from (for zero-pad columns there is none).
        plot(ax, full_d.t, full_d.drift, 'Color', [cols(v,:) 0.9], 'LineWidth', 1.6, ...
             'DisplayName', sprintf('%.2f s', delay_vals(v)));
        if strcmp(conds(c).type, 'delayed_start')
            xline(ax, delay_vals(v), '--', 'Color', cols(v,:), 'LineWidth', 1.4, ...
                  'HandleVisibility','off');
        end
    end
    xline(ax, 0, 'k-', 'LineWidth', 1.2, 'HandleVisibility','off');
    xlabel(ax, 'time from freeze onset (s)');
    if c == 1, ylabel(ax, 'drift  \beta\cdot sm(t)'); end
    title(ax, conds(c).name, 'FontWeight','normal');
    switch c
        case 1
            legend(ax, 'Location','northeast','Box','off','FontSize',12);
            note(ax, 'evidence DELETED on [0,L)');
        case 2, note(ax, 'real pre-onset evidence on [0,L)');
        case 3, note(ax, 'signal identical; dashed = accumulation start');
        case 4, note(ax, 'signal identical (curves overlap)');
    end
    apply_generic(ax, 'xlim', sig_xlim, 'ylim', sig_ylim, 'font_size', 15, 'line_width', 1.5);
end

% ---------- Row 2: example trajectories (same seeds everywhere) ----------
for c = 1:nC
    ax = nexttile(nC + c); hold(ax, 'on')
    for v = 1:numel(delay_vals)
        p = base; p.(conds(c).type) = delay_vals(v);
        if conds(c).use_pre, sig = sm_full; p.n_pre = n_pre;
        else,                sig = sm_post; p.n_pre = 0; end
        drift  = signal_to_drift(sig, p);          % already trimmed to t >= 0
        dstart = 0; if strcmp(conds(c).type, 'delayed_start'), dstart = delay_vals(v); end
        ndt_c  = 0; if strcmp(conds(c).type, 'ndt'),           ndt_c  = delay_vals(v); end
        n_skip = round(dstart / dt);
        drift_sim = drift(n_skip+1 : end);

        for s = 1:numel(seeds_traj)
            % x0 here is ABSOLUTE, matching what simulate_freeze_ddm passes to the
            % backends for row 3 (the wrapper takes the fraction instead).
            [k, traj] = sim_leaky_accumulator(drift_sim, theta, sigma, lambda, ...
                                              dt, seeds_traj(s), x0);
            t_abs = dstart + (1:numel(traj))' * dt;     % back to freeze-relative time
            plot(ax, t_abs, traj, '-', 'Color', [cols(v,:) 0.75], 'LineWidth', 1.4);
            if ~isnan(k)
                t_cross = dstart + k*dt;
                plot(ax, t_cross, theta, 'o', 'MarkerSize', 6, 'MarkerFaceColor', cols(v,:), ...
                     'MarkerEdgeColor','none');
                if ndt_c > 0   % ndt: path is identical, only the readout moves
                    plot(ax, [t_cross t_cross+ndt_c], [theta theta], ':', ...
                         'Color', cols(v,:), 'LineWidth', 1.6);
                    plot(ax, t_cross+ndt_c, theta, 'v', 'MarkerSize', 7, ...
                         'MarkerFaceColor', cols(v,:), 'MarkerEdgeColor','none');
                end
            end
        end
    end
    yline(ax, theta, 'k-', '\theta', 'LineWidth', 1.4, 'FontSize', 18);
    xlabel(ax, 'time from freeze onset (s)');
    if c == 1, ylabel(ax, 'DV'); end
    switch c
        case 1, note(ax, 'starved of evidence, then late');
        case 2, note(ax, 'same noise, slightly stale evidence');
        case 3, note(ax, 'path starts later');
        case 4, note(ax, 'paths identical; \nabla = reported RT');
    end
    apply_generic(ax, 'xlim', [0 xmax], 'ylim', traj_ylim, 'font_size', 15, 'line_width', 1.5);
end

% ---------- Row 3: full RT distributions ----------
% bins aligned to the dt grid (RTs are exact multiples of dt) to avoid aliasing
edges = -dt/2 : 3*dt : xmax;
for c = 1:nC
    ax = nexttile(2*nC + c); hold(ax, 'on')
    for v = 1:numel(delay_vals)
        p = base; p.(conds(c).type) = delay_vals(v);
        if conds(c).use_pre, sig = sm_full; p.n_pre = n_pre;
        else,                sig = sm_post; p.n_pre = 0; end
        rt = simulate_freeze_ddm(sig, p, N);
        RT{c, v} = rt;
        r = rt(~isnan(rt));
        histogram(ax, r, edges, 'Normalization','pdf', 'EdgeColor','none', ...
                  'FaceColor', cols(v,:), 'FaceAlpha', 0.35, ...
                  'DisplayName', sprintf('%.2f s  (med %.2f)', delay_vals(v), median(r)));
        xline(ax, median(r), '--', 'Color', cols(v,:), 'LineWidth', 1.6, ...
              'HandleVisibility','off');
    end
    xlabel(ax, 'reaction time (s)');
    if c == 1, ylabel(ax, 'pdf'); end
    legend(ax, 'Location','northeast','Box','off','FontSize',11);
    apply_generic(ax, 'xlim', [0 xmax], 'ylim', pdf_ylim, 'font_size', 15, 'line_width', 1.5);
end

%% ── Console summary: median shift per condition ─────────────────────────
fprintf('\nMedian RT (s) and shift vs the 0 s baseline:\n');
fprintf('%-34s', 'condition');
fprintf('%10.2f s', delay_vals); fprintf('\n');
M = nan(nC, numel(delay_vals));
for c = 1:nC
    M(c,:) = cellfun(@(r) median(r(~isnan(r))), RT(c,:));
    fprintf('%-34s', conds(c).name);
    fprintf('%12.3f', M(c,:)); fprintf('\n');
    fprintf('%-34s', '   shift vs baseline');
    fprintf('%12.3f', M(c,:) - M(c,1)); fprintf('\n');
end

% The headline: how much of the "sensory delay" effect is actually evidence loss?
fprintf('\nZero-pad vs real pre-onset stream (same sensory_delay):\n');
fprintf('%-12s %14s %14s %14s\n', 'delay', 'zero-pad', 'real pre-onset', 'due to padding');
for v = 1:numel(delay_vals)
    fprintf('%-12.2f %14.3f %14.3f %14.3f\n', delay_vals(v), M(1,v), M(2,v), M(1,v)-M(2,v));
end
fprintf(['\nAt delay 0 the two sensory columns are identical (internal control).\n' ...
         'Above it, the gap is the part of the "sensory delay" effect that is really\n' ...
         'EVIDENCE DELETION from zero-padding, not lag. Feeding the real pre-onset\n' ...
         'stream turns sensory_delay into a pure lag.\n']);

function o = ternary(c, a, b)
    if c, o = a; else, o = b; end
end
