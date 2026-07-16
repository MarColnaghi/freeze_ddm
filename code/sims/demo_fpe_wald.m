%% ════════════════════════════════════════════════════════════════════════
%  CRANK-NICOLSON FOKKER-PLANCK: a study demo
%  ------------------------------------------------------------------------
%  A readable MATLAB CN solver (fpe_cn.m -- read its header for the full
%  derivation) checked against two things:
%
%   (1) EXACT TRUTH. For a CONSTANT drift, no leak, and a single absorbing
%       barrier, the first-passage density is the Wald / inverse Gaussian
%       (wald_fpt.m). No grid, no sampling -- so any disagreement is solver
%       error, and we can watch it vanish as dx and dt shrink. This is a much
%       sharper instrument than comparing to Monte Carlo, which has its own
%       noise AND its own discretisation bias.
%
%   (2) A REAL FREEZE BOUT. The Wald only exists for constant drift; the whole
%       point of the PDE is a TIME-VARYING drive. So we then run the solver on
%       one real freeze's social-motion trace and cross-check three independent
%       implementations: this MATLAB solver, the C++ mex (leaky_pde_robust, via
%       ddm_pdf_from_trace), and the Monte-Carlo simulator.
%
%  Run headless:
%    matlab -batch "run('sims/demo_fpe_wald.m')"
% ════════════════════════════════════════════════════════════════════════
clearvars; clc;

this_dir = fileparts(mfilename('fullpath'));
addpath(fullfile(this_dir, 'simulators'));
addpath(fullfile(this_dir, 'simulators', 'cpp_mex_code'));
addpath(fullfile(fileparts(this_dir), 'fitting', 'instantaneous'));

% ── Regime ───────────────────────────────────────────────────────────────
dt       = 1/60;
dx       = 0.005;
sigma_sq = 1.0;   sigma = sqrt(sigma_sq);
x_min    = -10;   % far below: the Wald assumes NO lower boundary, so the
                  % solver's reflecting floor must be irrelevant
bound    = 1.5;
x0       = 0;
mu_const = 1.2;
T        = 631;
idx_trial = 4468;   % a real freeze with a clean, active motion trace

prm = struct('dt',dt, 'dx',dx, 'sigma_sq',sigma_sq, 'x0',x0, 'x_min',x_min, ...
             'bound',bound, 'lambda',0);

%% ── (1) CONSTANT DRIFT: solver vs the exact Wald ────────────────────────
fprintf('\n=== (1) CN solver vs exact Wald / inverse Gaussian ===\n');
drift_c = mu_const * ones(T,1);
[pdf_cn, surv_cn, t] = fpe_cn(drift_c, prm);

a = bound - x0;                              % barrier distance from the start
[pdf_w, cdf_w] = wald_fpt(t, a, mu_const, sigma);

% integrated comparison (CDF) is the fair one: the CN pdf is a per-step average,
% the Wald pdf is a point value, so they differ at O(dt) pointwise by construction.
cdf_cn = cumsum(pdf_cn) * dt;
ks     = max(abs(cdf_cn - cdf_w));
fprintf('  bound=%.2f  mu=%.2f  sigma=%.2f  ->  P(cross by %.1fs): CN %.4f | Wald %.4f\n', ...
    bound, mu_const, sigma, t(end), cdf_cn(end), cdf_w(end));
fprintf('  mass check: int(pdf)dt + surv = %.6f   (must be 1)\n', cdf_cn(end) + surv_cn);
fprintf('  KS(CN, Wald) = %.5f\n', ks);

%% ── (1b) convergence: does the solver error vanish with dx and dt? ──────
fprintf('\n=== (1b) convergence toward the exact solution ===\n');
fprintf('  refine dx (dt = 1/60):\n');
fprintf('  %10s %12s\n', 'dx', 'KS vs Wald');
for dxi = [0.04 0.02 0.01 0.005 0.0025]
    q = prm; q.dx = dxi;
    [pd, ~, tq] = fpe_cn(drift_c, q);
    [~, cw]     = wald_fpt(tq, a, mu_const, sigma);
    fprintf('  %10.4f %12.5f\n', dxi, max(abs(cumsum(pd)*dt - cw)));
end
fprintf('  refine dt (dx = 0.005):\n');
fprintf('  %10s %12s\n', 'dt', 'KS vs Wald');
for dti = [1/30 1/60 1/120 1/240]
    q = prm; q.dt = dti;
    Ti = round(T*dt/dti);
    [pd, ~, tq] = fpe_cn(mu_const*ones(Ti,1), q);
    [~, cw]     = wald_fpt(tq, a, mu_const, sigma);
    fprintf('  %10.5f %12.5f\n', dti, max(abs(cumsum(pd)*dti - cw)));
end

%% ── (2) REAL FREEZE BOUT: three implementations, time-varying drive ─────
fprintf('\n=== (2) real freeze bout (trial %d): CN vs mex vs Monte-Carlo ===\n', idx_trial);
paths      = path_generator('folder', 'sims/sanity_checks');
bouts      = importdata(fullfile(paths.dataset, 'bouts.mat'));
bouts_proc = data_parser_new(bouts, 'type','immobility', 'period','loom', ...
                             'window','le', 'nloom', 2:20);
bouts_proc.smp       = ones(height(bouts_proc), 1);
bouts_proc.intercept = ones(height(bouts_proc), 1);
sm_signal  = extract_sm_from_bouts(bouts_proc, 'type','onsets', ...
                             'output_type','mat', 'window', [0 630], 'norm_factor', 10);
sm = fillmissing(sm_signal(idx_trial,:).', 'previous');
sm = fillmissing(sm, 'next');

beta_gain = 3.0;
drift_tv  = beta_gain * sm;

% (a) this MATLAB solver
[pdf_tv, surv_tv] = fpe_cn(drift_tv, prm);
cdf_tv = cumsum(pdf_tv)*dt;

% (b) the C++ mex, through ddm_pdf_from_trace (same scheme, same conventions)
fixed = struct('dt',dt, 'dx',dx, 'sigma_sq',sigma_sq, 'x0',x0, 'x_min',x_min, ...
               'T_max',T*dt, 'T_trunc',0);
o       = ddm_pdf_from_trace([0, bound, 0], drift_tv, fixed);
cdf_mex = cumsum(o.pdf_raw(:))*dt;

% (c) Monte-Carlo, at the PLAIN bound. Two conventions must be reconciled first
% (see sims/test_sim_vs_fpe.m):
%   - drift: the solver reads the END-of-step drift, the simulator the START, so
%     shift the simulator's copy one frame earlier;
%   - discrete monitoring: the simulator only tests x >= theta at frame edges, so
%     it behaves like a CONTINUOUS process with the barrier RAISED by
%     0.5826*sigma*sqrt(dt). The correction therefore belongs on the PDE's bound,
%     not the simulator's -- putting it on the simulator doubles the gap.
N_mc = 500000;
p_mc = struct('dt',dt, 'theta',bound, 'sigma',sigma, ...
              'leak',0, 'x0',x0, 'gain',1, 'bias',0, 'backend','mex', 'seed',11);
rt   = simulate_freeze_ddm([drift_tv(2:end); drift_tv(end)], p_mc, N_mc);
k    = round(rt/dt); k = k(~isnan(k)); k = k(k <= T-1);
cdf_mc = cumsum(accumarray(k+1, 1, [T 1]) / N_mc);

prm_bgk = prm; prm_bgk.bound = bound + 0.5826*sigma*sqrt(dt);
pdf_bgk = fpe_cn(drift_tv, prm_bgk);
cdf_bgk = cumsum(pdf_bgk)*dt;

fprintf('  mass check (CN): int(pdf)dt + surv = %.6f\n', cdf_tv(end) + surv_tv);
fprintf('  KS(CN, mex)              = %.3e  <- same scheme, must be ~0\n', max(abs(cdf_tv - cdf_mex)));
fprintf('  KS(CN, MC) naive         = %.4f   <- discrete-monitoring bias\n', max(abs(cdf_tv - cdf_mc)));
fprintf('  KS(CN+BGK, MC) corrected = %.4f   <- should fall to ~MC noise (%.4f)\n', ...
    max(abs(cdf_bgk - cdf_mc)), 1/sqrt(N_mc));

%% ── Figure ───────────────────────────────────────────────────────────────
figure('Color','w','Position',[40 40 1500 820]);
tiledlayout(2, 3, 'TileSpacing','compact', 'Padding','compact');

ax = nexttile; hold(ax,'on')
plot(ax, t, pdf_w,  '-',  'LineWidth', 2.6, 'Color', [0.20 0.45 0.75], 'DisplayName','Wald (exact)');
plot(ax, t, pdf_cn, '--', 'LineWidth', 1.4, 'Color', [0.90 0.35 0.13], 'DisplayName','CN solver');
xlabel(ax,'first-passage time (s)'); ylabel(ax,'density');
title(ax,'constant drift: solver vs exact','FontWeight','normal');
legend(ax,'Location','northeast','Box','off','FontSize',11);
apply_generic(ax,'xlim',[0 5],'font_size',14,'line_width',1.3);

ax = nexttile; hold(ax,'on')
plot(ax, t, cdf_w,  '-',  'LineWidth', 2.6, 'Color', [0.20 0.45 0.75]);
plot(ax, t, cdf_cn, '--', 'LineWidth', 1.4, 'Color', [0.90 0.35 0.13]);
xlabel(ax,'first-passage time (s)'); ylabel(ax,'CDF');
subtitle(ax, sprintf('KS = %.5f', ks), 'FontAngle','italic','Color',[0.35 0.35 0.35],'FontSize',11);
title(ax,'CDF','FontWeight','normal');
apply_generic(ax,'xlim',[0 5],'ylim',[0 1],'font_size',14,'line_width',1.3);

ax = nexttile; hold(ax,'on')
plot(ax, t, pdf_cn - pdf_w, '-', 'LineWidth', 1.4, 'Color', [0.3 0.3 0.3]);
yline(ax, 0, 'k:');
xlabel(ax,'first-passage time (s)'); ylabel(ax,'CN - Wald (density)');
title(ax,'residual','FontWeight','normal');
apply_generic(ax,'xlim',[0 5],'font_size',14,'line_width',1.3);

ax = nexttile; hold(ax,'on')
plot(ax, t, drift_tv, '-', 'LineWidth', 1.2, 'Color', [0.2 0.5 0.3]);
xlabel(ax,'time from freeze onset (s)'); ylabel(ax,'drift  \beta\cdot sm(t)');
title(ax, sprintf('real freeze bout (trial %d)', idx_trial),'FontWeight','normal');
apply_generic(ax,'xlim',[0 10.5],'font_size',14,'line_width',1.3);

ax = nexttile; hold(ax,'on')
plot(ax, t, pdf_tv,        '-',  'LineWidth', 2.6, 'Color', [0.20 0.45 0.75], 'DisplayName','CN (MATLAB)');
plot(ax, t, o.pdf_raw(:),  '--', 'LineWidth', 1.4, 'Color', [0.90 0.35 0.13], 'DisplayName','mex (C++)');
xlabel(ax,'first-passage time (s)'); ylabel(ax,'density');
title(ax,'time-varying drive','FontWeight','normal');
legend(ax,'Location','northeast','Box','off','FontSize',11);
apply_generic(ax,'xlim',[0 10.5],'font_size',14,'line_width',1.3);

ax = nexttile; hold(ax,'on')
plot(ax, t, cdf_tv,  '-',  'LineWidth', 2.6, 'Color', [0.20 0.45 0.75], 'DisplayName','CN (plain bound)');
plot(ax, t, cdf_bgk, '-',  'LineWidth', 1.6, 'Color', [0.20 0.65 0.45], 'DisplayName','CN + BGK bound');
plot(ax, t, cdf_mc,  ':',  'LineWidth', 2.2, 'Color', [0.20 0.20 0.20], 'DisplayName','Monte-Carlo');
xlabel(ax,'first-passage time (s)'); ylabel(ax,'CDF');
title(ax,'CN vs Monte-Carlo','FontWeight','normal');
subtitle(ax, 'MC matches the BGK-corrected bound, not the plain one', ...
         'FontAngle','italic','Color',[0.35 0.35 0.35],'FontSize',11);
legend(ax,'Location','southeast','Box','off','FontSize',11);
apply_generic(ax,'xlim',[0 10.5],'ylim',[0 1],'font_size',14,'line_width',1.3);
