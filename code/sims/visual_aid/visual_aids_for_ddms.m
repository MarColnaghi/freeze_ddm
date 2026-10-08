%% ════════════════════════════════════════════════════════════════════════
%  VISUAL AID: what the drift gain does to a few DDM trajectories.
%  ------------------------------------------------------------------------
%  Three accumulator paths (three fixed seeds) drawn against the bound, then
%  redrawn at a lower beta -- the figure is built once and the line/marker data
%  are swapped in place, so the change in slope is the only thing that moves.
%
%  Only the paths that are actually drawn get simulated (3 per beta, not 2000),
%  they are integrated as one matrix instead of frame by frame, and the social-
%  motion pipeline -- ~11 s of loading and parsing -- is skipped unless the
%  drift really is driven by it. Same seeds, same figure, ~40x less work.
% ════════════════════════════════════════════════════════════════════════
clearvars; close all; clc;
paths = path_generator('folder', 'sims/visual_aid');

dt    = 1/60;                 % s / frame (native social-motion rate)
sigma = 1.0;                  % diffusion SD
theta = 3.1;                  % bound -- keep SMALL: at the extrema end (checks
                              % 4/5) a crossing must happen within ONE frame,
                              % where the noise SD is only sigma*sqrt(dt) ~ 0.13,
                              % so a large theta would make the extrema detector
                              % ~never fire (endpoint -> 100% censored).
leak  = 1;                    % 0 = perfect accumulator
x0    = 2.4;                    % initial position, ABSOLUTE units (same scale as
                              % theta, matching sim_leaky_accumulator) -- pass
                              % e.g. 0.3*theta for a bayes_fpe-style head start
betas = [1.3];            % drift gains, shown one after the other
idx_traj = [42 44 450];       % seeds of the paths on show
n = 631;                      % frames = width of the [0 630] social-motion window

%% ── Drive ───────────────────────────────────────────────────────────────
%  A constant drive is what the figure has always drawn (the old script loaded
%  the real signal and then overwrote it), so it is the default and costs
%  nothing. Set use_real_drift = true to drive the accumulator with a REAL
%  per-frame social-motion bout instead, loaded with the same pipeline as the
%  validation section at the bottom of this file.
use_real_drift = false;
sm = ones(n, 1);

if use_real_drift
    paths      = path_generator('folder', 'sims/visual_aid');
    idx_trial  = 950; % 2021; 5222; 11; 3394; 2020; 1918; 1001; 100 is great;
                      % 4548; 878; 917; 6666; 3945; 412; 7327; 3331 is fantastic
    bouts      = importdata(fullfile(paths.dataset, 'bouts.mat'));
    bouts_proc = data_parser_new(bouts, 'type','immobility', 'period','loom', ...
                                 'window','le', 'nloom', 2:20);
    sm_signal  = extract_sm_from_bouts(bouts_proc, 'type','onsets', ...
                                 'output_type','mat', 'window', [0 630], 'norm_factor', 10);
    sm = fillmissing(sm_signal(idx_trial, :).', 'previous');
    n  = numel(sm);
end

%% ── Figure ──────────────────────────────────────────────────────────────
fh = figure('Color','w','Position',[10 10 600 310]);
tl = tiledlayout(1, 1, 'TileSpacing', 'tight', 'Padding', 'tight');
hold on
dot_size = 260;
frames   = (0:n-1)';

for b = 1:numel(betas)

    [traj, last_frame] = accumulate(betas(b)*sm, theta, sigma, leak, dt, idx_traj, x0);

    if b == 1
        trj_hndl = plot(frames, traj, 'LineWidth', 1.5, 'Color', [0.75 0.75 0.75]);
        plot([-50 n + 50], [theta theta], 'LineWidth', 3.5, 'Color', 'k');
        plot([0 0], [theta 0], 'k--', 'LineWidth', 3.5);
        plot([0 0], [0 -3], 'k-', 'LineWidth', 3.5);

        scatter(0, x0, dot_size, 'filled', 'k')

        sc_hndl = gobjects(1, numel(idx_traj));
        for i = 1:numel(idx_traj)
            sc_hndl(i) = scatter(last_frame(i), theta, dot_size, 'filled', ...
                                 'MarkerFaceColor', [1 0.33 0.33]);
        end

        apply_generic(gca, 'xlim', [-25 n + 25], 'ylim', [-.75 theta + 0.3], 'no_y', true, 'no_x', true)
    else
        pause(1)
        for i = 1:numel(idx_traj)
            trj_hndl(i).YData = traj(:,i);
            sc_hndl(i).XData  = last_frame(i);
        end
    end
  %  exporter(fh, paths, sprintf('ddm_schematic_leak%g_x0%g_mu_%d.pdf', leak, x0, b))
end

%% ── The distribution those three paths came from ────────────────────────
%  The two conditions of the schematics: same leak, same drift, one WITHOUT a
%  head start and one starting at x0 = 2.4, just under the bound. With leak = 1
%  these are the distributions behind ddm_schematic_leak1_mu_1 and
%  ddm_schematic_leak1_x02.4_mu_1.
%
%  Defective density: normalised by ALL trials, so the area under a curve is the
%  fraction that crossed inside the window and the two are comparable. Log y,
%  because the early spike (crossings off the head start) is ~100x the shelf
%  (escapes from the well) and would otherwise flatten it.
x0_dist = [0 1.9];            % the two start points compared
b_dist  = 1;                  % which beta stage: 1 -> betas(1), the mu_1 figures
n_dist  = 20000;
leak = .7;
theta = 2.7;
bw      = 0.25;  edges = dt/2:dt*5:n*dt;  ctr = edges(1:end-1) + bw/2;
cols_d  = [0 0 0; 1 0.33 0.33];

fh2 = figure('Color','w','Position',[640 10 600 345]); hold on
for c = 1:numel(x0_dist)
    rt = rt_only(betas(b_dist)*sm, theta, sigma, leak, dt, x0_dist(c), n_dist);
    histogram(rt, edges, 'EdgeColor', 'none')
end
xlabel('Duration'); ylabel('density')
apply_generic(gca, 'xlims', [-0.2 n*dt + 0.2], 'ylim', [])
exporter(fh2, paths, sprintf('ddm_rtdist_leak%g_x0_%g_vs_%g.pdf', leak, x0_dist))

%% ════════════════════════════════════════════════════════════════════════
function [traj, last_frame] = accumulate(drift, theta, sigma, leak, dt, seeds, x0)
%ACCUMULATE  The paths of a handful of seeds, integrated as one matrix.
%   Same step, same noise stream and same crossing rule as sim_leaky_accumulator
%   -- one randn per frame, in frame order, under rng(seed) -- so column k is
%   exactly the path sim_leaky_accumulator(..., seeds(k)) returns. The loop runs
%   over the 3 seeds, not over the 631 frames.
%
%   x0 : initial position (optional, default 0), absolute units like theta.
%
%   traj       : n x numel(seeds), NaN after each path's crossing, with the
%                crossing sample snapped onto theta so the marker sits on it.
%   last_frame : crossing frame on the 0-based x axis, NaN for a path that never
%                crossed -- scatter() draws nothing at NaN, so an undecided path
%                simply gets no marker instead of a faked one at the last frame.

    n = numel(drift);
    m = numel(seeds);

    z = zeros(n, m);
    for k = 1:m
        rng(seeds(k));
        z(:,k) = randn(n, 1);
    end

    if nargin < 7 || isempty(x0), x0 = 0; end

    decay = 1 - leak*dt;
    inc   = drift(:)*dt + sigma*sqrt(dt)*z;
    if decay == 1
        traj = cumsum(inc, 1);                      % leak = 0: perfect integrator
    else
        traj = filter(1, [1 -decay], inc, [], 1);   % leaky: exp(-leak*t) memory
    end
    if x0 ~= 0                                      % a head start that leaks away
        traj = traj + x0 * decay.^((1:n)');         % too: x0 decays like the rest
    end

    % An absorbing bound leaves the path BEFORE the crossing untouched, which is
    % why the whole bundle can be integrated first and cut afterwards.
    [crossed, k_hit]    = max(traj >= theta, [], 1);
    crossed             = logical(crossed);
    last_frame          = nan(1, m);                % never crossed -> no decision,
    last_frame(crossed) = k_hit(crossed);           % so no crossing frame and no dot

    traj((1:n)' > last_frame) = NaN;                % NaN compares false: nothing cut
    traj(sub2ind(size(traj), k_hit(crossed), find(crossed))) = theta;

    last_frame = last_frame - 1;                    % 1-based frame -> 0-based axis
end

function rt = rt_only(drift, theta, sigma, leak, dt, x0, N)
%RT_ONLY  First-passage times only, for N trials.
%   Identical step and crossing rule to accumulate(), but the paths are thrown
%   away chunk by chunk instead of kept, so N can be large without the n x N
%   matrix. The seed is fixed here, so the two beta stages see the SAME noise
%   and differ only in drift. NaN = never crossed inside the window.

    n = numel(drift);  decay = 1 - leak*dt;  chunk = 5000;
    rt = nan(1, N);
    rng(0);
    for a = 1:chunk:N
        m = min(chunk, N - a + 1);
        x = filter(1, [1 -decay], drift(:)*dt + sigma*sqrt(dt)*randn(n, m), [], 1) ...
            + x0 * decay.^((1:n)');
        [cr, k] = max(x >= theta, [], 1);  cr = logical(cr);
        r = nan(1, m);  r(cr) = k(cr)*dt;
        rt(a:a+m-1) = r;
    end
end
