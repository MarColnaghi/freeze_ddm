%% ════════════════════════════════════════════════════════════════════════
%  ANATOMY OF A FIRST-PASSAGE DISTRIBUTION (single absorbing bound)
%  ------------------------------------------------------------------------
%  How the FPT density is actually computed, in six panels. The whole idea is
%  one line of arithmetic, and panels 4-6 are just three views of it:
%
%      S(t) = INT p(x,t) dx        <- survival: mass still un-crossed
%      f(t) = -dS/dt               <- the FPT density IS the mass lost per unit time
%
%  Why that works for a ONE-bound process: the lower edge only reflects (it is a
%  containment wall), so the ONLY way probability can leave the domain is by
%  being absorbed at the threshold. Every gram of mass that disappears is a
%  trial that just crossed. With TWO absorbing bounds this identity would lump
%  both together and you would have to track the flux at each bound separately.
%
%  Panel guide:
%    1. individual trajectories  -- what a trial is
%    2. the same thing as a DENSITY p(x,t) -- what the PDE tracks
%    3. snapshots: p spreads, drifts, and is PINNED TO ZERO at the bound
%    4. S(t): the mass draining away
%    5. f(t) = -dS/dt, against the exact Wald -- the payoff
%    6. F(t) + S(t) = 1: nothing is lost or invented
%
%  Run headless:
%    matlab -batch "run('sims/demo_fpe_anatomy.m')"
% ════════════════════════════════════════════════════════════════════════
clearvars; clc;

this_dir = fileparts(mfilename('fullpath'));
addpath(fullfile(this_dir, 'simulators'));
addpath(fullfile(this_dir, 'simulators', 'cpp_mex_code'));

% ── Regime: constant drift, no leak -> the exact Wald exists as ground truth
dt    = 1/60;   dx = 0.005;   sigma_sq = 1.0;  sigma = 1.0;
bound = 1.5;    x0 = 0;       mu = 1.2;        T  = 361;
x_min = -6;     % far below: keeps the Wald valid. We only PLOT the top of it.

prm = struct('dt',dt,'dx',dx,'sigma_sq',sigma_sq,'x0',x0,'x_min',x_min, ...
             'bound',bound,'lambda',0);
drift = mu*ones(T,1);

% Take the first 2 steps with backward Euler (Rannacher startup) so the density
% panels are readable. CN is A-stable but not L-stable -- its highest modes have
% amplification -> -1, so a delta start rings forever (the checkerboard). BE is
% L-stable and annihilates those modes; CN then carries on at O(dt^2). NB this
% changes only the TIME-STEPPING, not the model: the start is still a delta.
% (A Gaussian bump alone does NOT fix it -- it only cuts the ringing from
% min(p) = -185 to -26.)
prm_clean = prm; prm_clean.rannacher = 2;

[pdf, surv, t, x, P]  = fpe_cn(drift, prm_clean);   % P = full density p(x,t)
[~, ~, ~, ~, P_d]     = fpe_cn(drift, prm);         % plain CN: rings, for panel 3

S   = sum(P,2) * dx;                            % survival: INT p dx at each step
F   = cumsum(pdf) * dt;                         % CDF
[pdf_w, cdf_w] = wald_fpt(t, bound - x0, mu, sigma);

fprintf('mass check: F(end) + S(end) = %.6f   (must be 1)\n', F(end) + S(end));
fprintf('KS(CN, Wald) = %.5f\n', max(abs(F - cdf_w)));
fprintf('max |S(t) - (1 - F(t))| = %.2e   (survival and CDF are the same fact)\n', ...
    max(abs(S - (1 - F))));

snap_t = [0.15 0.4 0.8 1.6];                    % snapshot times (s)
snap_k = round(snap_t/dt) + 1;
cmapS  = [0.15 0.35 0.65; 0.20 0.60 0.55; 0.85 0.55 0.15; 0.75 0.25 0.25];
xl_plot = [-1.6 bound+0.15];

figure('Color','w','Position',[30 30 1620 860]);
tiledlayout(2,3,'TileSpacing','compact','Padding','compact');

%% 1. individual trajectories ---------------------------------------------
ax = nexttile; hold(ax,'on')
patch(ax, [0 t(end) t(end) 0], [bound bound xl_plot(2) xl_plot(2)], ...
      [0.95 0.9 0.9], 'EdgeColor','none');
for s = 1:14
    [k, traj] = sim_leaky_accumulator(drift, bound, sigma, 0, dt, s, x0);
    tt = (1:numel(traj))'*dt;
    plot(ax, tt, traj, '-', 'Color', [0.4 0.4 0.4 0.55], 'LineWidth', 0.8);
    if ~isnan(k)
        plot(ax, k*dt, bound, 'o', 'MarkerSize',5, 'MarkerFaceColor',[0.75 0.25 0.25], ...
             'MarkerEdgeColor','none');
    end
end
yline(ax, bound, 'k-', '\theta (absorbing)', 'LineWidth',1.6, 'FontSize',11);
yline(ax, x0, 'k:', 'x_0', 'LineWidth',1.2, 'FontSize',11);
xlabel(ax,'time (s)'); ylabel(ax,'accumulator x');
title(ax,'1. trials: each path stops when it first hits \theta','FontWeight','normal');
apply_generic(ax,'xlim',[0 2.5],'ylim',xl_plot,'font_size',13,'line_width',1.2);

%% 2. the ensemble as a density -------------------------------------------
ax = nexttile; hold(ax,'on')
imagesc(ax, t, x, P.');                          % p(x,t)
set(ax,'YDir','normal'); colormap(ax, flipud(gray(256)));
clim(ax, [0 0.9]);                               % the t=0 delta would blow the scale
for j = 1:numel(snap_k)
    xline(ax, t(snap_k(j)), '-', 'Color', cmapS(j,:), 'LineWidth', 1.6);
end
yline(ax, bound, 'r-', 'LineWidth', 2);
text(ax, 1.55, bound-0.16, 'mass is eaten here', 'Color',[0.75 0.1 0.1], ...
     'FontSize',11, 'FontAngle','italic');
xlabel(ax,'time (s)'); ylabel(ax,'x');
title(ax,'2. ...are one density p(x,t) the PDE evolves','FontWeight','normal');
apply_generic(ax,'xlim',[0 2.5],'ylim',xl_plot,'font_size',13,'line_width',1.2);

%% 3. snapshots ------------------------------------------------------------
ax = nexttile; hold(ax,'on')
% the same snapshot under PLAIN CN: high modes are undamped, so the delta start
% rings forever (min(p) reaches -185 here). Harmless for f(t) -- the lobes are
% symmetric, cancel in the mass integral, and sit far from theta -- but it is why
% p must not be clipped at 0. Rannacher startup removes it entirely.
plot(ax, x, P_d(snap_k(2),:), '-', 'LineWidth', 0.8, 'Color', [0.6 0.6 0.6], ...
     'DisplayName','plain CN (rings)');
for j = 1:numel(snap_k)
    plot(ax, x, P(snap_k(j),:), '-', 'LineWidth', 2, 'Color', cmapS(j,:), ...
         'DisplayName', sprintf('t = %.2f s', t(snap_k(j))));
end
xline(ax, bound, 'r-', 'LineWidth', 2, 'HandleVisibility','off');
text(ax, bound-0.62, 0.72, 'p(\theta,t) = 0', 'Color',[0.75 0.1 0.1], 'FontSize',12);
xlabel(ax,'x'); ylabel(ax,'p(x,t)');
title(ax,'3. it spreads, drifts, and is pinned to 0 at \theta','FontWeight','normal');
legend(ax,'Location','northwest','Box','off','FontSize',9);
apply_generic(ax,'xlim',xl_plot,'ylim',[0 0.85],'font_size',13,'line_width',1.2);

%% 4. survival -------------------------------------------------------------
ax = nexttile; hold(ax,'on')
plot(ax, t, S, '-', 'LineWidth', 2.6, 'Color', [0.2 0.2 0.2]);
for j = 1:numel(snap_k)
    plot(ax, t(snap_k(j)), S(snap_k(j)), 'o', 'MarkerSize',7, ...
         'MarkerFaceColor',cmapS(j,:), 'MarkerEdgeColor','none');
end
xlabel(ax,'time (s)'); ylabel(ax,'S(t) = \int p(x,t) dx');
title(ax,'4. mass drains out of the domain','FontWeight','normal');
subtitle(ax,'the only exit is \theta, so every loss is a crossing', ...
         'FontAngle','italic','Color',[0.35 0.35 0.35],'FontSize',11);
apply_generic(ax,'xlim',[0 2.5],'ylim',[0 1.02],'font_size',13,'line_width',1.2);

%% 5. the FPT density ------------------------------------------------------
ax = nexttile; hold(ax,'on')
plot(ax, t, pdf_w, '-', 'LineWidth', 3.0, 'Color', [0.75 0.80 0.90], ...
     'DisplayName','Wald (exact)');
plot(ax, t, pdf,   '-', 'LineWidth', 1.6, 'Color', [0.15 0.35 0.65], ...
     'DisplayName','-dS/dt from panel 4');
xlabel(ax,'first-passage time (s)'); ylabel(ax,'f(t)');
title(ax,'5. f(t) = -dS/dt  IS  the FPT density','FontWeight','normal');
subtitle(ax, sprintf('KS vs exact = %.5f', max(abs(F-cdf_w))), ...
         'FontAngle','italic','Color',[0.35 0.35 0.35],'FontSize',11);
legend(ax,'Location','northeast','Box','off','FontSize',11);
apply_generic(ax,'xlim',[0 2.5],'font_size',13,'line_width',1.2);

%% 6. bookkeeping ----------------------------------------------------------
ax = nexttile; hold(ax,'on')
area(ax, t, [F(:) S(:)], 'EdgeColor','none');
colormap(ax, [0.15 0.35 0.65; 0.85 0.85 0.85]);
plot(ax, t, F+S, 'k-', 'LineWidth', 1.6);
text(ax, 1.2, 0.35, 'F(t): already crossed', 'Color','w', 'FontSize',12, 'FontWeight','bold');
text(ax, 1.2, 0.82, 'S(t): still going', 'Color',[0.3 0.3 0.3], 'FontSize',12);
xlabel(ax,'time (s)'); ylabel(ax,'probability');
title(ax,'6. F(t) + S(t) = 1 at all times','FontWeight','normal');
subtitle(ax,'nothing lost, nothing invented','FontAngle','italic', ...
         'Color',[0.35 0.35 0.35],'FontSize',11);
apply_generic(ax,'xlim',[0 2.5],'ylim',[0 1.05],'font_size',13,'line_width',1.2);
