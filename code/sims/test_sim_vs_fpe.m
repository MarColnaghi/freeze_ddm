% TEST_SIM_VS_FPE  Validate the Monte-Carlo simulator against the FPE likelihood.
%
%   Every other check we have tests the SIGNAL FRONT-END (the bayes_fpe delay
%   port, the ReLU) or tests the accumulator against ITSELF (mex vs MATLAB).
%   This tests the accumulator DYNAMICS against an independent exact solution:
%   the Crank-Nicolson Fokker-Planck solver fitting/instantaneous/leaky_pde_robust.cpp,
%   wrapped by ddm_pdf_from_trace.m.
%
%     simulate_freeze_ddm  -> N sampled first-passage times   (stochastic)
%     ddm_pdf_from_trace   -> the first-passage DENSITY        (deterministic PDE)
%
%   They do NOT solve the same model, so a single-dt pass/fail would be
%   meaningless. The test is a CONVERGENCE argument instead:
%
%     a DISCRETISATION artefact vanishes as dt -> 0;  a MODEL difference does not.
%
%   Cases A-C differ only by discretisation, so their KS must shrink toward the
%   Monte-Carlo noise floor as dt shrinks. Case D differs by MODEL (see 3 below),
%   so its KS must PLATEAU. That contrast is the whole test.
%
%   The three known discrepancies:
%
%   1. DISCRETE-MONITORING BIAS (dominant). The PDE solves CONTINUOUS-time first
%      passage; the simulator is a discrete-time walk that only tests x >= theta
%      at frame boundaries, so it misses excursions that cross and return within
%      one dt and therefore crosses LATE. Classic Broadie-Glasserman-Kou: discrete
%      monitoring behaves like a continuous barrier raised by ~0.5826*sigma*sqrt(dt).
%      At dt=1/60, sigma=1 that is 0.075 -- only 2.8% of a bound of 2.7, but 21%
%      of a bound of 0.35. The naive gap scales as sqrt(dt) (ratio -> 1.414 in the
%      convergence table below), which is what identifies it.
%      => PRACTICAL CONSEQUENCE: fitting with the FPE and simulating with the MC
%         mixes two models. Expect the fitted bound to sit ~0.5826*sigma*sqrt(dt)
%         BELOW the bound you simulated with, unless you correct for it.
%
%   2. DRIFT-SAMPLING CONVENTION. leaky_pde_robust reads drifts[k] (0-based) for
%      the step ENDING at t=k*dt -- the drift at the END of the step. Our
%      sim_leaky_accumulator uses drift(k), the drift at the START. On a
%      time-varying drive that is a one-frame offset. Aligned here by prepending
%      one sample to the PDE's copy of the drift.
%
%   3. REFLECTING FLOOR (a real model difference). The PDE reflects at x_min; our
%      simulators have no floor at all (x runs to -inf). Irrelevant when x_min is
%      far from the action, decisive when it is not. Case D provokes it.
%
%   Verified separately: the PDE is CONVERGED IN SPACE at dx <= 0.005 (KS is
%   identical at dx = 0.005 / 0.0025 / 0.001), and the cell Peclet number stays
%   below ~1.5, so leaky_pde_robust's centred advection is not a source of error
%   here despite bayes_fpe preferring upwind for its leak branch.
%
%   Run headless:
%     matlab -batch "run('sims/test_sim_vs_fpe.m')"

this_dir = fileparts(mfilename('fullpath'));
addpath(fullfile(this_dir, 'simulators'));
addpath(fullfile(this_dir, 'simulators', 'cpp_mex_code'));
addpath(fullfile(fileparts(this_dir), 'fitting', 'instantaneous'));

sigma_sq = 1.0;
sigma    = sqrt(sigma_sq);
dx       = 0.005;          % converged (see header)
N        = 500000;
Ttot     = 6.0;
dt_list  = [1/60 1/120 1/240];
mc_noise = 1/sqrt(N);

% time-varying drive as a FUNCTION of t, so it can be resampled at each dt
drive = @(t) max(0.5 + 0.4*sin(2*pi*0.2*t) + 0.1*sin(2*pi*0.7*t), 0);

cases = struct( ...
  'name',    {'A constant drift',  'B time-varying',    'C leaky (lambda=5)', 'D near floor (x_min=-0.6)'}, ...
  'driftfun',{@(t) 2.0+0*t,        @(t) 4.2*drive(t),   @(t) 4.2*drive(t),    @(t) 0.35*drive(t)}, ...
  'bound',   {2.7,                 2.7,                 0.8,                  1.0}, ...
  'lambda',  {0,                   0,                   5,                    0}, ...
  'x_min',   {-7,                  -7,                  -7,                   -0.6}, ...
  'discret', {true,                true,                true,                 false});
%   discret = true  -> difference is discretisation only, KS must vanish as dt->0
%   discret = false -> difference is the MODEL (missing floor), KS must persist

fprintf('\nMC noise floor ~%.4f (N=%d), dx=%g\n\n', mc_noise, N, dx);
fprintf('%-28s %26s\n', '', 'KS (BGK-corrected + aligned)');
hdr = arrayfun(@(d) sprintf('dt=1/%d', round(1/d)), dt_list, 'UniformOutput', false);
fprintf('%-28s', 'case'); fprintf('%12s', hdr{:}); fprintf('%12s\n','naive@1/60');
KS = nan(numel(cases), numel(dt_list));
for ci = 1:numel(cases)
    C = cases(ci);
    fprintf('%-28s', C.name);
    for di = 1:numel(dt_list)
        KS(ci,di) = run_pair(C, dt_list(di), dx, sigma_sq, sigma, Ttot, N, true);
        fprintf('%12.4f', KS(ci,di));
    end
    ks_naive = run_pair(C, dt_list(1), dx, sigma_sq, sigma, Ttot, N, false);
    fprintf('%12.4f\n', ks_naive);
end

%% ---- convergence rate: sqrt(dt) identifies discrete monitoring -------------
fprintf('\nNaive gap vs dt (case A, no corrections) -- expect ratio -> sqrt(2)=1.414:\n');
fprintf('%10s %10s %10s\n','dt','KS','ratio');
prev = NaN;
for d = [1/60 1/120 1/240 1/480]
    ks = run_pair(cases(1), d, dx, sigma_sq, sigma, Ttot, N, false);
    if isnan(prev), fprintf('%10.5f %10.4f %10s\n', d, ks, '-');
    else,           fprintf('%10.5f %10.4f %10.3f\n', d, ks, prev/ks); end
    prev = ks;
end

%% ---- verdicts -------------------------------------------------------------
fprintf('\n');
tol_conv = 0.015;    % what "converged to the MC noise floor" means at dt=1/240
for ci = 1:numel(cases)
    C = cases(ci); ks_fine = KS(ci,end); ks_coarse = KS(ci,1);
    if C.discret
        ok = ks_fine < tol_conv && ks_fine < 0.6*ks_coarse;
        report(sprintf('%s: KS vanishes as dt->0 (%.4f -> %.4f)', ...
               C.name, ks_coarse, ks_fine), ok);
    else
        ok = ks_fine > 0.05 && ks_fine > 0.6*ks_coarse;
        report(sprintf('%s: KS PERSISTS as dt->0 (%.4f -> %.4f) = model difference', ...
               C.name, ks_coarse, ks_fine), ok);
    end
end
fprintf(['\nA-C converging to the MC noise floor validates the accumulator dynamics\n' ...
         '(drift, leak, bound, start point) against an independent exact solution.\n' ...
         'D refusing to converge is the point: no amount of dt refinement fixes a\n' ...
         'missing reflecting floor, because it is a different model.\n']);

%% ---- figure ---------------------------------------------------------------
nP = numel(cases);
figure('Color','w','Position',[30 30 1750 780]);
tiledlayout(2, nP, 'TileSpacing','compact', 'Padding','compact');
for ci = 1:nP
    [~, R] = run_pair(cases(ci), dt_list(1), dx, sigma_sq, sigma, Ttot, N, true);
    ax = nexttile(ci); hold(ax,'on')
    plot(ax, R.t, R.pmf_fpe, '-',  'LineWidth', 2.4, 'Color', [0.20 0.45 0.75], ...
         'DisplayName','FPE density');
    plot(ax, R.t, R.pmf_mc,  '--', 'LineWidth', 1.3, 'Color', [0.90 0.35 0.13], ...
         'DisplayName','Monte-Carlo');
    xlabel(ax,'first-passage time (s)'); if ci==1, ylabel(ax,'P(cross in frame)'); end
    title(ax, cases(ci).name, 'FontWeight','normal');
    if ci==1, legend(ax,'Location','northeast','Box','off','FontSize',11); end
    apply_generic(ax, 'xlim', [0 4], 'font_size', 14, 'line_width', 1.3);

    ax = nexttile(nP+ci); hold(ax,'on')
    plot(ax, dt_list, KS(ci,:), 'o-', 'LineWidth', 2, 'Color', [0.2 0.2 0.2], ...
         'MarkerFaceColor',[0.2 0.2 0.2]);
    yline(ax, mc_noise, ':', 'MC noise', 'LineWidth', 1.4, 'FontSize', 11);
    set(ax,'XScale','log','YScale','log');
    xlabel(ax,'dt (s)'); if ci==1, ylabel(ax,'KS (corrected)'); end
    if cases(ci).discret, s='discretisation -> vanishes'; else, s='MODEL diff -> persists'; end
    subtitle(ax, s, 'FontAngle','italic','Color',[0.35 0.35 0.35],'FontSize',11);
    apply_generic(ax, 'font_size', 14, 'line_width', 1.3);
end

%% ---- helpers --------------------------------------------------------------
function [ks, R] = run_pair(C, dt, dx, sigma_sq, sigma, Ttot, N, correct)
% correct=false : same drift, same bound for both (naive)
% correct=true  : align the drift convention AND apply the BGK bound correction
    T     = round(Ttot/dt);
    t     = (0:T-1)' * dt;
    drift = C.driftfun(t);
    x0    = C.x_min + round((0 - C.x_min)/dx)*dx;   % snap to the PDE grid

    drift_fpe = drift;
    bound_fpe = C.bound;
    if correct
        % the mex reads drifts[k] (END of step) where we use drift(k) (START)
        drift_fpe = [drift(1); drift(1:end-1)];
        % discrete monitoring == continuous barrier raised by 0.5826*sigma*sqrt(dt)
        bound_fpe = C.bound + 0.5826*sigma*sqrt(dt);
    end

    fixed = struct('dt',dt,'dx',dx,'sigma_sq',sigma_sq,'x0',x0, ...
                   'x_min',C.x_min,'T_max',T*dt,'T_trunc',0);
    o = ddm_pdf_from_trace([C.lambda, bound_fpe, 0], drift_fpe, fixed);

    p  = struct('dt',dt,'theta',C.bound,'sigma',sigma,'leak',C.lambda,'x0',x0, ...
                'gain',1,'bias',0,'backend','mex','seed',11);
    rt = simulate_freeze_ddm(drift, p, N);

    k = round(rt/dt); k = k(~isnan(k)); k = k(k <= T-1);
    R.pmf_mc  = accumarray(k+1, 1, [T 1]) / N;
    R.pmf_fpe = o.pdf_raw(:) * dt;
    R.cdf_mc  = cumsum(R.pmf_mc);
    R.cdf_fpe = cumsum(R.pmf_fpe);
    R.t       = o.t(:);
    ks = max(abs(R.cdf_mc - R.cdf_fpe));
end

function report(name, cond)
    if cond, tag='PASS'; else, tag='FAIL'; end
    fprintf('  [%s] %s\n', tag, name);
end
