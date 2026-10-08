%% ════════════════════════════════════════════════════════════════════════
%  WHAT THE LEAK DOES TO THE RT DISTRIBUTION
%  ------------------------------------------------------------------------
%  The same four settings as the trajectory schematics (ddm_schematic_leak*),
%  but with 40k paths instead of 3, so the shape behind those few lines is
%  visible. Every condition splits into TWO components:
%
%    1. an early spike  -- paths that cross straight off the head start at x0,
%                          before the leak has had time to forget it;
%    2. a flat tail     -- once the head start is gone the accumulator sits in
%                          its well at mu/leak, and the only way out is a noise
%                          excursion. Escapes have no preferred moment, so the
%                          hazard goes FLAT and the density decays exponentially
%                          (a straight line on the log axis of panel 1).
%
%  Panels, and why all three are here: the density alone is misleading under
%  this much censoring, because most trials never cross at all and so are not
%  in it. Panel 2 shows that directly -- where each curve ends is the fraction
%  still undecided when the window closes. Note it is still FALLING there: these
%  trials are right-censored by the 10.5 s window, not proven never to decide.
%  Panel 3 is the mechanism.
%
%  Densities are DEFECTIVE: normalised by the total number of trials, not by
%  the number that crossed, so the area under each curve is P(cross) and the
%  four curves are directly comparable in height.
%
%  Run headless:
%    matlab -batch "run('sims/visual_aid/ddm_rt_distributions.m')"
% ════════════════════════════════════════════════════════════════════════
clearvars; close all; clc;
paths = path_generator('folder', 'sims/visual_aid');

dt    = 1/60;                 % s / frame
sigma = 1.0;                  % diffusion SD
theta = 3.1;                  % bound
x0    = 2.4;                  % start point -- close to the bound, as in the schematics
n     = 631;                  % frames = the 10.5 s observation window
N     = 40000;                % trials per condition
chunk = 8000;                 % paths integrated at once (keeps memory flat)
seed  = 11;

%  leak, beta, and how to draw it. Colour = leak, line = drift gain.
col   = cmapper();
c_lo  = hex2rgb(col.processes.long);        % slower forgetting
c_hi  = hex2rgb(col.processes.short);       % faster forgetting
conds = { ...
    0.7, 1.3, c_lo, '-' ; ...
    1.0, 1.3, c_hi, '-' ; ...
    0.7, 0.5, c_lo, '--'; ...
    1.0, 0.5, c_hi, '--'};

bw    = 0.25;                               % s, histogram bin width
edges = 0:bw:n*dt;
ctr   = edges(1:end-1) + bw/2;

fh = figure('Color','w','Position',[20 20 1500 430]);
tl = tiledlayout(1, 3, 'TileSpacing','compact', 'Padding','compact');
ax = gobjects(1,3);
for a = 1:3, ax(a) = nexttile; hold(ax(a),'on'); end

lbl = cell(size(conds,1),1);
h1  = gobjects(size(conds,1),1);

for c = 1:size(conds,1)

    [leak, beta, colr, ls] = conds{c,:};
    rt = simulate_rts(beta, theta, sigma, leak, dt, n, x0, N, chunk, seed);

    crossed = ~isnan(rt);
    cnt     = histcounts(rt(crossed), edges);

    dens = cnt / (N * bw);                              % defective density
    surv = 1 - cumsum(cnt)/N;                           % P(RT > t), censored included
    haz  = cnt ./ max(N*[1 surv(1:end-1)] * bw, eps);   % crossings per survivor per s

    keep = surv > 0.01;                                 % stop where survivors run out
    h1(c) = plot(ax(1), ctr(dens>0), dens(dens>0), ls, 'Color', colr, 'LineWidth', 2.2);
            plot(ax(2), ctr, surv,  ls, 'Color', colr, 'LineWidth', 2.2);
            plot(ax(3), ctr(keep), haz(keep), ls, 'Color', colr, 'LineWidth', 2.2);

    lbl{c} = sprintf('\\lambda %.1f, \\beta %.1f  (%.0f%% by 10.5 s)', ...
                     leak, beta, 100*mean(crossed));
    fprintf(['leak %.1f  beta %.1f | ceiling %.2f | %5.1f%% by 10.5 s | ' ...
             'late hazard %.4f /s | median RT of those that DID cross %.2f s\n'], ...
            leak, beta, beta/leak, 100*mean(crossed), mean(haz(ctr>6 & keep)), ...
            median(rt(crossed)));
end

set(ax(1), 'YScale','log');  set(ax(3), 'YScale','log');
title(ax(1), 'RT density (defective)',        'FontWeight','normal','FontSize',15);
title(ax(2), 'survival: not yet decided',      'FontWeight','normal','FontSize',15);
title(ax(3), 'hazard: chance of deciding now','FontWeight','normal','FontSize',15);
ylabel(ax(1), 'density'); ylabel(ax(2), 'P(no decision yet)'); ylabel(ax(3), 'hazard (s^{-1})');
xlabel(tl, 'time since onset (s)', 'FontSize', 15);

legend(ax(1), h1, lbl, 'Location','northeast', 'Box','off', 'FontSize', 11);
ylim(ax(2), [0 1]);
for a = 1:3
    apply_generic(ax(a), 'xlims', [0 n*dt], 'font_size', 13, 'line_width', 1.2);
end

exporter(fh, paths, sprintf('ddm_rt_distributions_x0%g.pdf', x0))

%% ════════════════════════════════════════════════════════════════════════
function rt = simulate_rts(beta, theta, sigma, leak, dt, n, x0, N, chunk, seed)
%SIMULATE_RTS  First-passage times of N leaky accumulators, in blocks.
%   Same step and crossing rule as accumulate() in visual_aids_for_ddms.m; the
%   paths themselves are thrown away, only the crossing frame is kept, so the
%   memory cost is one chunk at a time rather than N x n.
%   rt is NaN for a trial that never crossed inside the window (right-censored).

    decay = 1 - leak*dt;
    rt    = nan(1, N);
    rng(seed);
    for a = 1:chunk:N
        b = min(a + chunk - 1, N);
        z = randn(n, b - a + 1);
        x = filter(1, [1 -decay], beta*dt + sigma*sqrt(dt)*z, [], 1) ...
            + x0 * decay.^((1:n)');
        [crossed, k] = max(x >= theta, [], 1);
        crossed = logical(crossed);
        r = nan(1, b - a + 1);
        r(crossed) = k(crossed) * dt;
        rt(a:b) = r;
    end
end
