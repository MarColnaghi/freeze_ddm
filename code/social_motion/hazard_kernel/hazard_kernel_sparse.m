clearvars; clc;
% ════════════════════════════════════════════════════════════════════════
% Sparse-prior hazard kernel with an EXPONENTIAL basis — companion to
% hazard_kernel_sm.m. Estimates the social-motion → un-freeze hazard kernel
% β(τ) three ways and horse-races them, following Mineault, Barthelmé & Pack
% (2009), "Improved classification images with sparse priors in a smooth basis":
%
%   (A) Smoothness / L2  — raised-cosine basis + P-spline (2nd-diff) penalty.
%                          This is the current hazard_kernel_sm.m method
%                          (the paper's Table-1 "smoothness" Gaussian prior).
%   (B) Sparse   / L1    — a BANK OF DECAYING EXPONENTIALS exp(-τ/τ_k).
%                          Under an L1 prior the fit becomes TIMESCALE
%                          SELECTION: the retained τ_k IS the leak estimate,
%                          testing integration (slow τ) vs transient (fast τ).
%   (C) Sparse   / L1    — flexible raised-cosine bumps (NEUTRAL comparator).
%
% Why (C): an exponential basis embeds the leaky-integrator hypothesis
% (paper's Criterion 2), so it must be adjudicated against a basis that can
% represent arbitrary smooth shapes. Cross-validated deviance decides:
%   (B) ≈ (A)/(C) at much lower effective df ⇒ kernel is a few exponentials
%       ⇒ integration with an identifiable leak = the selected τ_k.
%   (C) ≫ (B)                                 ⇒ non-exponential structure
%       exponentials miss.
% Significance of the kernel SHAPE is still judged by the time-shuffle null
% in hazard_kernel_sm.m — a sparse point estimate is not a significance test.
% ════════════════════════════════════════════════════════════════════════

fps = 60;

% ── Kernel window (relative to the un-freeze decision; 0 = escape) ───────
kernel_past_s   = 4.0;   % seconds of causal history
kernel_future_s = 0;     % seconds of acausal control (should come out ≈ 0)

% ── Discretisation / basis / fit settings ────────────────────────────────
dt_frames   = 6;      % hazard sampled every dt_frames (6 = 0.1 s → 10 Hz)
nb_kernel   = 12;     % raised-cosine bumps over CAUSAL lags (bases A, C)
nb_acausal  = 6;      % bumps over the ACAUSAL lags (all bases)
bump_width  = 1.25;   % raised-cosine half-width in units of centre spacing
nb_baseline = 6;      % raised-cosine bumps for the baseline hazard vs time-in-freeze
n_exp       = 10;      % number of exponential atoms in the bank (basis B)
tau_min_s   = 0.05;    % fastest exp time-constant (s)
tau_max_s   = 8.0;     % slowest exp time-constant (s); > kernel_past_s probes slow integration
                       % (atoms past the window are only weakly identified — see note in report)
cv_folds    = 10;      % grouped (by bout) cross-validation folds
save_out    = true;   % save results .mat + figure to disk
mask_preonset = true; % true → zero out lags reaching before freeze onset (no pre-onset drive)
grid_anchor = 'offset';% sampling-grid alignment: 'entry' = forward from onset (short bin on the
                      % break); 'offset' = backward from the freeze break (regular bin on the
                      % break, short bin at entry)
rng(120);

% ── Load + filter (same pipeline as hazard_kernel_sm.m) ───────────────────
threshold_imm = 2; threshold_mob = 2; threshold_pc = 4;
id_code = sprintf('imm%d_mob%d_pc%d', threshold_imm, threshold_mob, threshold_pc);
paths        = path_generator('folder', 'social_motion/hazard_kernel', 'bouts_id', id_code, 'imfirst', false);
bouts        = importdata(fullfile(paths.dataset, 'bouts.mat'));
motion_cache = importdata(fullfile(paths.cache_path, 'motion_cache.mat'));

thresholds = define_thresholds('le_window', struct('le_window_sl', [5 55], 'le_window_fl', [5 55]));
bouts      = bouts_formatting(bouts, thresholds);

trunc_point  = 30;
bl = data_parser_new(bouts, 'type', 'immobility', 'period', 'loom', 'window', 'le', ...
        'nloom', 2:20, 'min_dur', trunc_point);

contact_threshold = 70;   % px
bl = impose_contact_threshold(bl, 'threshold', contact_threshold, 'type', 'onlyfreeze');

bl.dur_frames  = bl.ending_time;
bl.bout_id     = (1:height(bl))';
bl = bl(bl.dur_frames > trunc_point, :);
entry_fr = trunc_point;   % delayed entry (left truncation)

fprintf('Bouts: %d total, %d censored by collision (%.1f%%)\n', ...
    height(bl), sum(bl.is_censored), 100*mean(bl.is_censored));
fprintf('Settings: mask_preonset=%d  grid_anchor=%s\n', mask_preonset, grid_anchor);

% ── Lag grid ──────────────────────────────────────────────────────────────
back_fr =  round(kernel_past_s   * fps);   % positive lags = past
fwd_fr  =  round(kernel_future_s * fps);   % negative lags = future (control)
lag_fr  = (-fwd_fr : back_fr)';            % τ in frames; s(τ)=sm(t−τ)
nLag    = numel(lag_fr);
caus_mask = lag_fr >= 0;                   % causal (past): τ ≥ 0
t_rel   = -lag_fr / fps;                   % time rel. to escape (s): <0 past, >0 future

% ── Global z-scoring of sm ────────────────────────────────────────────────
flies = unique(bl.fly);
nn = 0; s1 = 0; s2 = 0;
for i = 1:numel(flies)
    v = motion_cache(flies(i)); v = v(:); v = v(~isnan(v));
    nn = nn + numel(v); s1 = s1 + sum(v); s2 = s2 + sum(v.^2);
end
sm_mu = s1/nn;  sm_sd = sqrt(s2/nn - sm_mu^2);

% ── Build long (person-period) design: raw history S + outcome + nuisances ─
max_rows = sum(ceil(bl.dur_frames / dt_frames)) + height(bl);
S    = nan(max_rows, nLag);        % raw z-scored sm history (projected per-basis below)
tinf = nan(max_rows, 1);           % time in freeze (s) — baseline hazard
yev  = zeros(max_rows, 1);         % 1 = genuine un-freeze event this interval
gid  = nan(max_rows, 1);           % bout id (for grouped CV)
mov  = nan(max_rows, 1);           % moving_flies condition
slo  = nan(max_rows, 1);           % loom speed

nB = height(bl);
fprintf('Building person-period design (%d bouts)...\n', nB);
t0 = tic;
r = 0;
for b = 1:nB
    sm  = motion_cache(bl.fly(b)); sm = sm(:);
    L   = numel(sm);
    on  = bl.onsets(b);
    dur = bl.dur_frames(b);

    if dur < entry_fr, continue; end
    % Sampling grid on [entry_fr, dur], stepping dt_frames. Anchored either at
    % entry (forward from onset; irregular bin lands on the break) or at the
    % break (backward from offset; regular bin on the break, irregular bin at entry).
    switch grid_anchor
        case 'entry',  base = entry_fr;
        case 'offset', base = dur;
        otherwise, error('grid_anchor must be ''entry'' or ''offset''.');
    end
    kk   = ceil((1-base)/dt_frames) : floor((dur-base)/dt_frames);
    allg = base + kk(:)*dt_frames;
    grid = allg(allg >= entry_fr);
    if strcmp(grid_anchor,'entry')
        if isempty(grid) || grid(end) < dur,     grid = [grid; dur];      end %#ok<AGROW>
    else
        if isempty(grid) || grid(1)  > entry_fr, grid = [entry_fr; grid]; end %#ok<AGROW>
    end

    for gi = 1:numel(grid)
        f     = grid(gi);
        t_abs = on + f - 1;
        if (t_abs - back_fr) < 1 || (t_abs + fwd_fr) > L, continue; end

        s = (sm(t_abs - lag_fr) - sm_mu) / sm_sd;   % z-scored history, length nLag
        if any(isnan(s)), continue; end
        if mask_preonset, s(lag_fr > f-1) = 0; end  % drop pre-onset drive (0 = mean, neutral)

        r = r + 1;
        S(r,:)  = s';
        tinf(r) = f / fps;
        gid(r)  = bl.bout_id(b);
        mov(r)  = bl.moving_flies(b);
        slo(r)  = bl.sloom(b);

        if gi == numel(grid) && ~bl.is_censored(b)   % event on final interval of uncensored bout
            yev(r) = 1;
        end
    end
    if mod(b, ceil(nB/10)) == 0
        fprintf('  design %3.0f%%  (%d/%d bouts, %d rows so far) [%.0fs]\n', 100*b/nB, b, nB, r, toc(t0));
    end
end

S=S(1:r,:); tinf=tinf(1:r); yev=yev(1:r); gid=gid(1:r); mov=mov(1:r); slo=slo(1:r);
fprintf('Design: %d at-risk intervals, %d un-freeze events (hazard %.3f) [%.0fs]\n', ...
    r, sum(yev), mean(yev), toc(t0));

% ── Baseline-hazard basis on time-in-freeze (equal-mass knots; drop one) ──
ut = tiedrank(tinf) / numel(tinf);
Bb = raised_cosine_basis(ut, nb_baseline);
Bb = Bb(:, 2:end);
nb_base_use = size(Bb, 2);

% ── Nuisance dummies ──────────────────────────────────────────────────────
[~,~,im] = unique(mov);  Dmov = dummyvar(im);  Dmov = Dmov(:,2:end);
[~,~,is] = unique(slo);  Dslo = dummyvar(is);  Dslo = Dslo(:,2:end);
Nuis     = [Bb, Dmov, Dslo];             % all unpenalised (besides intercept)
n_nuis   = size(Nuis, 2);

% ── Grouped CV folds (bout-level), shared across all models ───────────────
fold    = mod(randperm(max(gid)), cv_folds) + 1;
rowfold = fold(gid);

% ── Build the three kernel bases on lag_fr ────────────────────────────────
tau_grid_s = logspace(log10(tau_min_s), log10(tau_max_s), n_exp);   % fast → slow time-constants
[Bexp, tau_grid, is_exp] = exp_basis(lag_fr, tau_grid_s, fps, nb_acausal, bump_width);
Bflex = raised_cosine_basis(lag_fr, nb_kernel + nb_acausal, bump_width);
Bflex = Bflex ./ max(vecnorm(Bflex), eps);          % unit-norm cols (scale-fair L1)
Bsm   = raised_cosine_basis(lag_fr, nb_kernel + nb_acausal, bump_width);   % L2 baseline basis

fprintf('\nBases: exp bank = %d atoms (+%d acausal), flex = %d bumps, smooth = %d bumps\n', ...
    n_exp, nb_acausal, size(Bflex,2), size(Bsm,2));

% ════════════════════════════════════════════════════════════════════════
% (A) Smoothness / L2 on the raised-cosine basis  — current method
% ════════════════════════════════════════════════════════════════════════
Xk_sm   = S * Bsm;
XA      = [ones(r,1), Xk_sm, Nuis];
kerA    = 1 + (1:size(Bsm,2));
Dk      = diffmat(size(Bsm,2), 2);
PkA     = zeros(size(XA,2)); PkA(kerA,kerA) = Dk' * Dk;

lam_L2  = logspace(-1, 6, 12);
fprintf('\n[A] L2 smoothness: CV over %d lambdas x %d folds... [%.0fs]\n', numel(lam_L2), cv_folds, toc(t0));
[lambda_A, cvllA] = select_lambda_l2(XA, yev, PkA, lam_L2, rowfold, cv_folds);
[bA, VA] = penalized_logit(XA, yev, lambda_A * PkA);

kernel_A    = Bsm * bA(kerA);
Cov_kA      = VA(kerA, kerA);
kernel_A_se = sqrt(sum((Bsm * Cov_kA) .* Bsm, 2));
[devA, dfA, aicA] = model_metrics_l2(XA, yev, bA, VA);
cv_devA = -2 * cvllA;

% ════════════════════════════════════════════════════════════════════════
% (0) Nuisance-only baseline (kernel ABSENT) — the reference every kernel
%     model must beat. If the L1 fits collapse onto this, the sparse prior
%     is telling you the kernel adds no cross-validated predictive value.
% ════════════════════════════════════════════════════════════════════════
fprintf('\n[0] nuisance-only baseline... [%.0fs]\n', toc(t0));
Xnull = [ones(r,1), Nuis];
PtinyN = 1e-6 * eye(size(Xnull,2));
[bN, VN] = penalized_logit(Xnull, yev, PtinyN);
[devN, dfN, aicN] = model_metrics_l2(Xnull, yev, bN, VN);
cv_devN = -2 * pen_grouped_cv(Xnull, yev, PtinyN, rowfold, cv_folds);

% Selection rule for the L1 fits. 1-SE zeroes the (weak) kernel here, so we
% use argmax CV to let the kernel show and let exp vs flex actually diverge.
sel_rule = 'argmax';

% ════════════════════════════════════════════════════════════════════════
% (B) Sparse / L1 on the EXPONENTIAL bank
% ════════════════════════════════════════════════════════════════════════
fprintf('\n[B] L1 exponential bank: CV + full lambda path... [%.0fs]\n', toc(t0));
[bB, kernel_B, kernel_B_se, sel_tau, lambda_B, cvllB, devB, dfB, aicB, kerB] = ...
    fit_l1(S, Bexp, Nuis, yev, rowfold, cv_folds, is_exp, tau_grid, sel_rule, 'B: L1 exp');
cv_devB = -2 * cvllB;

% ════════════════════════════════════════════════════════════════════════
% (C) Sparse / L1 on the flexible bumps (neutral comparator)
% ════════════════════════════════════════════════════════════════════════
fprintf('\n[C] L1 flexible bumps: CV + full lambda path... [%.0fs]\n', toc(t0));
[bC, kernel_C, kernel_C_se, ~, lambda_C, cvllC, devC, dfC, aicC, kerC] = ...
    fit_l1(S, Bflex, Nuis, yev, rowfold, cv_folds, false(1,size(Bflex,2)), [], sel_rule, 'C: L1 flex');
cv_devC = -2 * cvllC;
fprintf('\nAll fits complete [%.0fs].\n', toc(t0));

% ════════════════════════════════════════════════════════════════════════
% Report
% ════════════════════════════════════════════════════════════════════════
fprintf('\n── Model comparison (lower AIC / lower CV deviance = better) ──\n');
fprintf('%-26s %8s %6s %10s %14s\n', 'model', 'dev', 'df', 'AIC', 'CV dev/obs');
fprintf('%-26s %8.1f %6.1f %10.1f %14.4f\n', '(0) nuisance-only',   devN, dfN, aicN, cv_devN);
fprintf('%-26s %8.1f %6.1f %10.1f %14.4f\n', '(A) L2 smoothness',   devA, dfA, aicA, cv_devA);
fprintf('%-26s %8.1f %6.1f %10.1f %14.4f\n', '(B) L1 exponential',  devB, dfB, aicB, cv_devB);
fprintf('%-26s %8.1f %6.1f %10.1f %14.4f\n', '(C) L1 flexible',     devC, dfC, aicC, cv_devC);
fprintf('(kernel adds CV value iff A/B/C beat row 0 by more than a few units)\n');
fprintf('\nλ:  A(L2)=%.3g   B(L1 exp)=%.3g   C(L1 flex)=%.3g\n', lambda_A, lambda_B, lambda_C);

if isempty(sel_tau)
    fprintf('(B) selected NO exponential atoms (kernel ≈ 0).\n');
else
    fprintf('(B) selected time-constants τ_k (s): %s\n', ...
        strjoin(arrayfun(@(t) sprintf('%.2f', t), sel_tau, 'uni', 0), ', '));
    fprintf('    → leak read-out: slow τ ⇒ integration, fast τ ⇒ transient.\n');
end

% ── Figure: overlay the three recovered kernels ──────────────────────────
fh = figure('Color','w','Position',[100 100 820 520]); hold on
fill([t_rel; flipud(t_rel)], [kernel_B + 1.96*kernel_B_se; flipud(kernel_B - 1.96*kernel_B_se)], ...
    [0.85 0.4 0.2], 'FaceAlpha', 0.15, 'EdgeColor','none')
pA = plot(t_rel, kernel_A, 'Color',[0.2 0.5 0.8], 'LineWidth', 2);
pB = plot(t_rel, kernel_B, 'Color',[0.85 0.4 0.2], 'LineWidth', 2);
pC = plot(t_rel, kernel_C, '--', 'Color',[0.2 0.5 0.8], 'LineWidth', 2);
yline(0,'k:'); xline(0,'--k')
% shade the acausal control region
yl = ylim;
patch([0 kernel_future_s kernel_future_s 0], [yl(1) yl(1) yl(2) yl(2)], ...
    [0.9 0.9 0.9], 'FaceAlpha', 0.35, 'EdgeColor','none')
uistack(findobj(gca,'Type','patch'),'bottom')
text(kernel_future_s/2, yl(2)*0.9, 'acausal control', 'HorizontalAlignment','center', 'FontSize',9)
xlabel('Time relative to un-freeze (s)')
ylabel('Kernel weight  \beta(\tau)')
title('Social-motion → un-freeze hazard kernel: L2 smooth vs L1 exponential vs L1 flexible')
legend([pA pB pC], {sprintf('(A) L2 smooth (df=%.0f)',dfA), ...
                    sprintf('(B) L1 exp (df=%.0f)',dfB), ...
                    sprintf('(C) L1 flex (df=%.0f)',dfC)}, 'Location','best', 'box','off')
apply_generic(gca)

% ── Persist ───────────────────────────────────────────────────────────────
if save_out
    stamp  = datestr(now, 'yyyymmdd_HHMMSS'); %#ok<TNOW1,DATST>
    outdir = paths.dataset;
    try
        params = struct('kernel_past_s',kernel_past_s, 'kernel_future_s',kernel_future_s, ...
            'dt_frames',dt_frames, 'nb_kernel',nb_kernel, 'nb_acausal',nb_acausal, ...
            'bump_width',bump_width, 'nb_baseline',nb_baseline, 'n_exp',n_exp, ...
            'cv_folds',cv_folds, 'contact_threshold',contact_threshold, 'trunc_point',trunc_point, ...
            'mask_preonset',mask_preonset, 'grid_anchor',grid_anchor);
        save(fullfile(outdir, ['hazard_kernel_sparse_' stamp '.mat']), ...
            'kernel_A','kernel_A_se','kernel_B','kernel_B_se','kernel_C','kernel_C_se', ...
            't_rel','lag_fr','caus_mask','fps','tau_grid','is_exp','sel_tau', ...
            'lambda_A','lambda_B','lambda_C', ...
            'devA','dfA','aicA','cv_devA','devB','dfB','aicB','cv_devB', ...
            'devC','dfC','aicC','cv_devC','bA','bB','bC','params', '-v7.3');
        exportgraphics(fh, fullfile(outdir, ['hazard_kernel_sparse_' stamp '.png']), 'Resolution',200);
        fprintf('\nSaved results + figure to %s (stamp %s)\n', outdir, stamp);
    catch ME
        warning('save_out failed: %s', ME.message);
    end
end

% ════════════════════════════════════════════════════════════════════════
% Local functions
% ════════════════════════════════════════════════════════════════════════
function [B, tau_grid, is_exp] = exp_basis(lag_fr, tau_grid_s, fps, nb_acausal, bump_width)
% Causal dictionary of decaying exponentials exp(-τ/τ_k) on τ≥0, plus a small
% raised-cosine block on the acausal (future) lags as a leakage control. Columns are
% L2-normalised so the L1 penalty is scale-fair across atoms (Mineault et al. 2009, Eq.11).
    lag_s = lag_fr(:) / fps;
    caus  = lag_s >= 0;
    ne    = numel(tau_grid_s);
    Be    = zeros(numel(lag_s), ne);
    for k = 1:ne
        col = zeros(numel(lag_s),1);
        col(caus) = exp(-lag_s(caus) / tau_grid_s(k));
        Be(:,k) = col / max(norm(col), eps);
    end
    Ba = zeros(numel(lag_s), nb_acausal);
    if any(~caus) && nb_acausal > 0
        Ba(~caus,:) = raised_cosine_basis(lag_fr(~caus), nb_acausal, bump_width);
        Ba = Ba ./ max(vecnorm(Ba), eps);
    end
    B        = [Be, Ba];
    tau_grid = tau_grid_s;
    is_exp   = [true(1,ne), false(1,nb_acausal)];
end

function [b, kernel, kernel_se, sel_tau, lambda, cvll, dev, df, aic, ker_cols] = ...
        fit_l1(S, Bk, Nuis, y, rowfold, K, is_exp_cols, tau_grid, sel_rule, name)
% Fit an L1-penalised logistic kernel in basis Bk, print the full grouped-CV λ-path
% (with the number of surviving kernel atoms at each λ), select λ by 'argmax' (best CV)
% or '1se' (largest λ within 1 SE), reconstruct the kernel with a debiased-refit CI, and
% return goodness-of-fit metrics.
    if nargin < 9 || isempty(sel_rule), sel_rule = 'argmax'; end
    if nargin < 10, name = 'L1'; end
    r        = size(S,1);
    Xk       = S * Bk;
    X        = [ones(r,1), Xk, Nuis];
    nb_k     = size(Bk,2);
    ker_cols = 1 + (1:nb_k);
    penmask  = false(1, size(X,2));  penmask(ker_cols) = true;

    % λ path: λ_max exactly zeros every penalised coef (KKT at the nuisance-only fit —
    % the working correlation collapses to |Xⱼᵀ(y−μ₀)| on the first IRLS step).
    b_nuis   = penalized_logit(X(:, ~penmask), y, 0);
    mu0      = 1 ./ (1 + exp(-min(max(X(:,~penmask)*b_nuis,-30),30)));
    lam_max  = max(abs(X(:,penmask)' * (y - mu0)));
    lam_grid = lam_max * logspace(0, -3, 12);        % descending: all-zero → nearly unpenalised

    % grouped-CV log-likelihood with warm starts along the path (per fold)
    llf = nan(K, numel(lam_grid));
    for k = 1:K
        fprintf('    [%s] CV fold %d/%d...\n', name, k, K);
        te = rowfold == k; tr = ~te;
        if ~any(te) || ~any(tr), continue; end
        bw = [];
        for il = 1:numel(lam_grid)
            bw = penalized_logit_l1(X(tr,:), y(tr), lam_grid(il), penmask, bw);
            eta = min(max(X(te,:)*bw, -30), 30);
            mu  = min(max(1 ./ (1 + exp(-eta)), 1e-9), 1-1e-9);
            yy  = y(te);
            llf(k,il) = mean(yy.*log(mu) + (1-yy).*log(1-mu));
        end
    end
    cvll_lam = mean(llf, 1, 'omitnan');
    cvse_lam = std(llf, 0, 1, 'omitnan') ./ sqrt(sum(~isnan(llf),1));

    % full-data path: fit at every λ (warm-started) and count surviving kernel atoms
    fprintf('    [%s] full-data lambda path...\n', name);
    natom = nan(1,numel(lam_grid));  Bfull = cell(1,numel(lam_grid));  bw = [];
    for il = 1:numel(lam_grid)
        bw = penalized_logit_l1(X, y, lam_grid(il), penmask, bw);
        Bfull{il} = bw;
        natom(il) = sum(abs(bw(ker_cols)) > 1e-8);
    end

    [~, iarg] = max(cvll_lam);                        % best CV
    thr   = cvll_lam(iarg) - cvse_lam(iarg);
    i1se  = find(cvll_lam >= thr, 1, 'first');        % largest λ (descending path) within 1 SE
    if strcmpi(sel_rule, '1se'), isel = i1se; else, isel = iarg; end
    lambda = lam_grid(isel);
    cvll   = cvll_lam(isel);
    b      = Bfull{isel};

    % print the path so the kernel's marginal predictive value is visible
    fprintf('\n── %s  λ-path  (* = argmax CV, s = 1-SE pick; sel=%s) ──\n', name, sel_rule);
    fprintf('%3s %11s %13s %7s\n', '#', 'lambda', 'CVll/obs', 'natom');
    for il = 1:numel(lam_grid)
        mk = '  '; if il==iarg, mk(1)='*'; end; if il==i1se, mk(2)='s'; end
        fprintf('%3d %11.3g %13.5f %7d %s\n', il, lam_grid(il), cvll_lam(il), natom(il), mk);
    end

    bk     = b(ker_cols);
    kernel = Bk * bk;

    % goodness-of-fit: deviance, df = #nonzero coefs (paper Eq. B6), AIC = dev + 2 df
    eta = min(max(X*b, -30), 30);  mu = min(max(1./(1+exp(-eta)),1e-9),1-1e-9);
    dev = -2 * sum(y.*log(mu) + (1-y).*log(1-mu));
    df  = sum(abs(b) > 1e-8);
    aic = dev + 2*df;

    % debiased-refit CI: unpenalised logistic on [nuisance | selected kernel atoms]
    active = ker_cols(abs(bk) > 1e-8) - 1;           % selected kernel-column indices (into Bk)
    if isempty(active)
        kernel_se = zeros(size(kernel));  sel_tau = [];
    else
        Xdb = [ones(r,1), S*Bk(:,active), Nuis];
        [~, Vdb] = penalized_logit(Xdb, y, 1e-6*eye(size(Xdb,2)));
        kcol     = 1 + (1:numel(active));
        Cov_k    = Vdb(kcol, kcol);
        Bact     = Bk(:, active);
        kernel_se = sqrt(sum((Bact * Cov_k) .* Bact, 2));
        if ~isempty(tau_grid) && any(is_exp_cols)
            sel_tau = tau_grid(active(is_exp_cols(active)));
        else
            sel_tau = [];
        end
    end
end

function [lambda, cvll] = select_lambda_l2(X, y, Pk, lam_grid, rowfold, K)
% Pick the L2 (P-spline) penalty strength by grouped CV, 1-SE rule (largest λ within 1 SE).
    cvll_lam = nan(size(lam_grid));  cvse_lam = nan(size(lam_grid));
    for il = 1:numel(lam_grid)
        [cvll_lam(il), llf] = pen_grouped_cv(X, y, lam_grid(il)*Pk, rowfold, K);
        cvse_lam(il) = std(llf, 'omitnan') / sqrt(sum(~isnan(llf)));
        fprintf('    [A] lambda %2d/%d = %9.3g  CVll/obs=%.5f\n', il, numel(lam_grid), lam_grid(il), cvll_lam(il));
    end
    [~, ibest] = max(cvll_lam);
    thr   = cvll_lam(ibest) - cvse_lam(ibest);
    i1se  = find(cvll_lam >= thr, 1, 'last');        % lam_grid ascending → largest λ within 1 SE
    lambda = lam_grid(i1se);
    cvll   = cvll_lam(i1se);
end

function [dev, df, aic] = model_metrics_l2(X, y, b, V)
% Deviance, EFFECTIVE df = tr((XᵀWX+P)⁻¹ XᵀWX) (paper Eq. B7), AIC = dev + 2 df.
    eta = min(max(X*b, -30), 30);  mu = min(max(1./(1+exp(-eta)),1e-9),1-1e-9);
    dev = -2 * sum(y.*log(mu) + (1-y).*log(1-mu));
    w   = max(mu.*(1-mu), 1e-9);
    XWX = X' * (X .* w);
    df  = trace(V * XWX);                            % V = inv(XWX + P)
    aic = dev + 2*df;
end

function B = raised_cosine_basis(x, nb, width_mult, anchor)
% Linearly-spaced raised-cosine bumps tiling the range of x (see hazard_kernel_sm.m).
    if nargin < 3 || isempty(width_mult), width_mult = 1; end
    if nargin < 4, anchor = []; end
    x = x(:);
    lo = min(x); hi = max(x);
    c  = linspace(lo, hi, nb);
    sp = (hi - lo) / (nb - 1);
    if ~isempty(anchor)
        [~, jn] = min(abs(c - anchor));
        c = c + (anchor - c(jn));
    end
    w  = width_mult * sp;
    B  = zeros(numel(x), nb);
    for j = 1:nb
        d = (x - c(j)) / w;
        in = abs(d) <= 1;
        B(in, j) = 0.5 * (1 + cos(pi * d(in)));
    end
end

function D = diffmat(n, ord)
% k-th order difference operator: (n-ord) x n. diffmat(n,2) penalises curvature.
    D = eye(n);
    for k = 1:ord, D = diff(D); end
end

function [b, V] = penalized_logit(X, y, P)
% Penalised (L2) logistic regression by IRLS / Newton; V = inv(XᵀWX + P).
    [~, p] = size(X);
    b   = zeros(p, 1);
    rdg = 1e-8 * eye(p);
    for it = 1:100 %#ok<NASGU>
        eta = min(max(X*b, -30), 30);
        mu  = 1 ./ (1 + exp(-eta));
        w   = max(mu .* (1 - mu), 1e-9);
        H   = X' * (X .* w) + P + rdg;
        g   = X' * (y - mu) - P * b;
        step = H \ g;
        b   = b + step;
        if max(abs(step)) < 1e-8, break; end
    end
    V = inv(H);
end

function [b, support] = penalized_logit_l1(X, y, lambda, penmask, b0)
% MAP logistic regression with an L1 penalty of strength lambda on the columns flagged by
% penmask (logical; 1 = penalised). Outer IRLS + inner cyclic coordinate descent with
% soft-thresholding. Unpenalised columns are updated by plain weighted least squares.
%
% Covariance-update (glmnet-style) form: per IRLS step we form the p×p weighted Gram
% A = XᵀWX and c = XᵀWz ONCE (one BLAS call, O(n·p²)), then run coordinate descent
% entirely in p-space by maintaining rc = c − A·b and updating it in O(p) per coordinate.
% Algebraically identical to residual-space CD (same fixed point), but each coordinate
% update is O(p) instead of O(n) — the win is largest when the columns are highly
% correlated (e.g. the exponential bank), where CD needs many sweeps to converge.
    [~, p] = size(X);
    if nargin < 5 || isempty(b0), b = zeros(p,1); else, b = b0(:); end
    penmask = logical(penmask(:));
    for outer = 1:100 %#ok<NASGU>
        eta = min(max(X*b, -30), 30);
        mu  = 1 ./ (1 + exp(-eta));
        w   = max(mu .* (1 - mu), 1e-6);
        z   = eta + (y - mu) ./ w;                 % IRLS working response
        Xw  = X .* w;                              % n×p (weighted design)
        A   = Xw' * X;                             % p×p weighted Gram  (== XᵀWX)
        c   = Xw' * z;                             % p×1  (== XᵀWz)
        d   = max(diag(A), eps);                   % Σ w x_j^2 per column
        b_old = b;
        rc  = c - A*b;                             % residual correlation  c − A·b
        for inner = 1:1000 %#ok<NASGU>
            bmax = 0;
            for j = 1:p
                rho = rc(j) + d(j)*b(j);           % = Σ w x_j (z − Σ_{k≠j} x_k b_k)
                if penmask(j)
                    bj = sign(rho) * max(abs(rho) - lambda, 0) / d(j);
                else
                    bj = rho / d(j);
                end
                db = bj - b(j);
                if db ~= 0, rc = rc - A(:,j)*db; b(j) = bj; bmax = max(bmax, abs(db)); end
            end
            if bmax < 1e-7, break; end
        end
        if max(abs(b - b_old)) < 1e-7, break; end
    end
    support = find(abs(b) > 1e-8);
end

function [ll, ll_folds] = pen_grouped_cv(X, y, P, rowfold, K)
% Held-out Bernoulli log-likelihood/obs of the penalised (L2) model, folds grouped by bout.
    ll_folds = nan(K, 1);
    for k = 1:K
        te = rowfold == k; tr = ~te;
        if ~any(te) || ~any(tr), continue; end
        b   = penalized_logit(X(tr,:), y(tr), P);
        eta = min(max(X(te,:)*b, -30), 30);
        mu  = min(max(1 ./ (1 + exp(-eta)), 1e-9), 1-1e-9);
        yy  = y(te);
        ll_folds(k) = mean(yy.*log(mu) + (1-yy).*log(1-mu));
    end
    ll = mean(ll_folds, 'omitnan');
end
