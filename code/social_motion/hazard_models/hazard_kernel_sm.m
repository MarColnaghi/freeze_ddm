clearvars; clc;
% ════════════════════════════════════════════════════════════════════════
% Time-resolved hazard model of freeze termination as a function of the
% social-motion signal the focal fly experiences.
%
% Instead of collapsing sm(t) to one scalar per bout, we model the
% moment-to-moment probability of un-freezing (the discrete-time hazard) as
% a function of the RECENT HISTORY of social motion, via a temporal kernel.
%
% Stage 1 — Temporal kernel (agnostic):
%   logit P(unfreeze at t | frozen) = baseline(time_in_freeze)
%                                      + Σ_τ β(τ)·sm(t−τ) + controls
%   The SHAPE of β(τ) distinguishes the computations:
%     • broad, sustained, decaranying slowly  → integration / accumulation
%       (exponential decay ⇒ leak time-constant of a leaky integrator)
%     • sharp, concentrated near τ≈0        → instantaneous / transient
%   A small acausal (future) segment is included as a control: it should be
%   flat ≈ 0. Weight there flags reverse causation or autocorrelation.
%
% Stage 2 — Integration vs extrema (explicit model comparison):
%   Extrema detection is nonlinear, so we also build scalar features over
%   the causal window — mean (integration) vs max (peak/extrema) vs slope
%   (transient) — and compare hazard models by AIC and bout-grouped
%   cross-validated log-likelihood. Whichever feature carries the
%   predictive weight names the computation.
%
% Censoring: a collision truncates the bout (true duration is LONGER), so
% those bouts are right-censored — they contribute at-risk intervals but
% never an un-freeze event.
% ════════════════════════════════════════════════════════════════════════

fps = 60;

% ── Kernel window (relative to the un-freeze decision; 0 = escape) ───────
% Causal history = social motion BEFORE the decision. Acausal = control.
kernel_past_s   = 3.0;   % seconds of causal history (set the integration horizon)
kernel_future_s = 1;   % seconds of acausal control (should come out ≈ 0)

% ── Discretisation / basis / fit settings ───────────────────────────────
dt_frames   = 6;     % hazard sampled every dt_frames (6 = 0.1 s → 10 Hz); cuts autocorrelation
nb_kernel   = 12;     % raised-cosine bumps over the CAUSAL (past) lags (lowered 5→4: fewer DOF ⇒ less ringing)
nb_acausal  = 4;     % bumps over the ACAUSAL (future) lags
bump_width  = 1.25;   % kernel bump half-width in units of center spacing (>1 ⇒ more overlap ⇒ smoother band)
% Kernel-basis + shuffle mode:
%   'full'  → single CONTINUOUS basis over the whole window + full-window shuffle.
%             Clean, continuous kernel; flat leakage-free null; no τ=0 seam. (recommended)
%   'split' → block-diagonal causal/acausal bases (no bump crosses τ=0) + per-half
%             shuffle. Lets you assess causal vs acausal ordering separately, but
%             introduces an artificial discontinuity / edge artifact at the break point.
% Both use nb_kernel+nb_acausal total bumps, so model complexity matches across modes.

shuffle_mode = 'split';
use_penalty = true;  % true → P-spline roughness penalty on the kernel (λ by 1-SE rule);
                     % false → plain unpenalised logistic kernel (original behaviour)
anchor_zero = false;  % 'full' mode: force one bump center exactly on τ=0, so a causal
                     % effect at escape isn't pushed to the nearest off-zero center
nb_baseline = 6;     % raised-cosine bumps for the baseline hazard vs time-in-freeze
cv_folds    = 5;     % grouped (by bout) cross-validation folds
n_shuffle   = 200;  % within-window time-shuffles (overnight: 1000+ → tight null & p-floor 1/(n+1))
save_out    = true;  % save results .mat + figures to disk (overnight runs)
% Per-bout-mean centering (Mundlak within/between decomposition): subtract each
% bout's own mean sm so the KERNEL reflects WITHIN-bout fluctuations only, and add
% the per-bout mean as a separate covariate that absorbs/reports the BETWEEN-bout
% level effect. Tests whether a causal kernel survives removing between-bout level
% (high-motion bouts being short) — i.e. is the timing effect genuinely within-bout.
center_per_bout = true;
rng(1);

% ── Load + filter (same pipeline as check_fd_distributions.m) ────────────
threshold_imm = 2; threshold_mob = 2; threshold_pc = 4;
id_code = sprintf('imm%d_mob%d_pc%d', threshold_imm, threshold_mob, threshold_pc);
paths        = path_generator('folder', 'social_motion/proximity', 'bouts_id', id_code, 'imfirst', false);
bouts        = importdata(fullfile(paths.dataset, 'bouts.mat'));
motion_cache = importdata(fullfile(paths.cache_path, 'motion_cache.mat'));

thresholds = define_thresholds('le_window', struct('le_window_sl', [5 55], 'le_window_fl', [5 55]));
bouts      = bouts_formatting(bouts, thresholds);

trunc_point  = 18;
bl = data_parser_new(bouts, 'type', 'immobility', 'period', 'loom', 'window', 'le', ...
        'nloom', 2:20, 'min_dur', trunc_point);

contact_threshold = 75;   % px — conservative, above the empirical contamination distance
bl = impose_contact_threshold(bl, 'threshold', contact_threshold, 'type', 'onlyfreeze');

% ending_time = observed (possibly censored) freeze length, in frames
bl.dur_frames  = bl.ending_time;
bl.bout_id     = (1:height(bl))';
bl = bl(bl.dur_frames > trunc_point, :);

% Delayed-entry time: the min-duration scoring rule is LEFT TRUNCATION.
% Bouts enter the risk set only once they have survived to entry_fr (0.5 s);
% see the loop below for why counting earlier intervals would bias the result.
entry_fr = trunc_point;

fprintf('Bouts: %d total, %d censored by collision (%.1f%%)\n', ...
    height(bl), sum(bl.is_censored), 100*mean(bl.is_censored));

% ── Lag grid + bases ─────────────────────────────────────────────────────
back_fr =  round(kernel_past_s   * fps);   % positive lags = past
fwd_fr  =  round(kernel_future_s * fps);   % negative lags = future (control)
lag_fr  = (-fwd_fr : back_fr)';            % τ in frames; s(τ)=sm(t−τ)
nLag    = numel(lag_fr);
caus_mask = lag_fr >= 0;                   % causal (past): τ ≥ 0
caus_t    = lag_fr(caus_mask) / fps;       % seconds into the past

anchor_lag = []; if anchor_zero, anchor_lag = 0; end   % knot on τ=0 (lag 0) when requested
switch shuffle_mode
    case 'full'
        % Single continuous basis over the whole window — smooth kernel, no τ=0 seam.
        Bk = raised_cosine_basis(lag_fr, nb_kernel + nb_acausal, bump_width, anchor_lag);
    case 'split'
        % Block-diagonal: causal bumps only on causal lags, acausal only on acausal
        % lags. No bump crosses τ=0, so a within-half shuffle is leakage-free — but
        % the split allows a discontinuity/edge artifact at the break point.
        Bk = zeros(nLag, nb_kernel + nb_acausal);
        Bk(caus_mask,  1:nb_kernel)              = raised_cosine_basis(lag_fr(caus_mask),  nb_kernel,  bump_width);
        Bk(~caus_mask, nb_kernel+(1:nb_acausal)) = raised_cosine_basis(lag_fr(~caus_mask), nb_acausal, bump_width);
    otherwise
        error('shuffle_mode must be ''full'' or ''split''');
end
nb_k = size(Bk, 2);

% ── Global z-scoring of sm (one scale; preserves kernel shape) ───────────
flies = unique(bl.fly);
n = 0; s1 = 0; s2 = 0;
for i = 1:numel(flies)
    v = motion_cache(flies(i)); v = v(:);
    v = v(~isnan(v));
    n = n + numel(v); s1 = s1 + sum(v); s2 = s2 + sum(v.^2);
end
sm_mu = s1/n;  sm_sd = sqrt(s2/n - sm_mu^2);

% ── Build long (person-period) design ────────────────────────────────────
max_rows = sum(ceil(bl.dur_frames / dt_frames)) + height(bl);
Xk    = nan(max_rows, nb_k);        % kernel-basis projections of sm history
S     = nan(max_rows, nLag);        % raw z-scored sm history (kept for the shuffle control)
tinf  = nan(max_rows, 1);           % time in freeze (s) — baseline hazard
yev   = zeros(max_rows, 1);         % 1 = genuine un-freeze event this interval
gid   = nan(max_rows, 1);           % bout id (for grouped CV)
mov   = nan(max_rows, 1);           % moving_flies condition
slo   = nan(max_rows, 1);           % loom speed
fmean = nan(max_rows, 1);           % Stage-2: mean sm over causal window
fmax  = nan(max_rows, 1);           % Stage-2: max  sm over causal window
fslp  = nan(max_rows, 1);           % Stage-2: slope of sm over causal window
bmean = nan(max_rows, 1);           % per-bout mean sm (Mundlak between-bout covariate)

% precompute slope regressor (linear fit of s_caus vs caus_t)
A_slope = [caus_t - mean(caus_t)];
A_slope = A_slope / (A_slope' * A_slope);   % gives slope = A_slope' * s_caus

r = 0;
for b = 1:height(bl)
    sm  = motion_cache(bl.fly(b)); sm = sm(:);
    L   = numel(sm);
    on  = bl.onsets(b);
    dur = bl.dur_frames(b);

    % Left truncation / delayed entry: start the risk set at entry_fr, not 0.
    % Counting intervals in (0, entry_fr) would create an immortal-time zone —
    % survivors present, but all sub-0.5 s terminations deleted upstream — which
    % drives the baseline hazard ≈ 0 there and, because high social motion tends
    % to end freezes fast, selectively drops high-sm events and distorts the
    % early kernel. So the bout contributes intervals only from entry_fr onward.
    if dur < entry_fr, continue; end
    grid = (entry_fr:dt_frames:dur)';
    if isempty(grid) || grid(end) < dur, grid = [grid; dur]; end %#ok<AGROW>

    % per-bout mean sm (z-scored), over the frozen frames — the between-bout level
    seg_idx = on : min(on + dur - 1, L);
    mbar    = (mean(sm(seg_idx), 'omitnan') - sm_mu) / sm_sd;
    if isnan(mbar), mbar = 0; end           % degenerate bout → no centering shift

    for gi = 1:numel(grid)
        f     = grid(gi);
        t_abs = on + f - 1;

        % require the full kernel window to be in-bounds
        if (t_abs - back_fr) < 1 || (t_abs + fwd_fr) > L, continue; end

        s = (sm(t_abs - lag_fr) - sm_mu) / sm_sd;   % z-scored history, length nLag
        if any(isnan(s)), continue; end
        if center_per_bout, s = s - mbar; end       % within-bout transform (Mundlak)

        r = r + 1;
        Xk(r,:) = (s' * Bk);                 % project onto kernel basis
        S(r,:)  = s';                        % keep raw history for the shuffle control
        tinf(r) = f / fps;
        gid(r)  = bl.bout_id(b);
        mov(r)  = bl.moving_flies(b);
        slo(r)  = bl.sloom(b);
        bmean(r)= mbar;                      % between-bout level covariate (const within bout)

        sc       = s(caus_mask);
        fmean(r) = mean(sc);
        fmax(r)  = max(sc);
        fslp(r)  = A_slope' * (sc - mean(sc));

        % event only on the final interval of an UNCENSORED bout
        if gi == numel(grid) && ~bl.is_censored(b)
            yev(r) = 1;
        end
    end
end

% trim
Xk=Xk(1:r,:); S=S(1:r,:); tinf=tinf(1:r); yev=yev(1:r); gid=gid(1:r);
mov=mov(1:r); slo=slo(1:r); fmean=fmean(1:r); fmax=fmax(1:r); fslp=fslp(1:r);
bmean=bmean(1:r);
fprintf('Design: %d at-risk intervals, %d un-freeze events (hazard %.3f)\n', ...
    r, sum(yev), mean(yev));

% Baseline-hazard basis on time-in-freeze.
% (1) Place knots on the empirical-CDF (rank) scale of tinf so every bump
%     covers an equal-mass slice of durations. Tiling raw time leaves the
%     rare long-duration bumps nearly empty → near-zero columns → rank
%     deficiency, worst inside CV folds. (2) Drop one bump: raised cosines
%     sum to ≈1 and would otherwise duplicate the intercept fitglm adds.
ut = tiedrank(tinf) / numel(tinf);          % empirical CDF in (0,1]
Bb = raised_cosine_basis(ut, nb_baseline);
Bb = Bb(:, 2:end);
nb_base_use = size(Bb, 2);

% ── Grouped CV folds (bout-level) — reused for λ-selection, Stage 2, shuffle ─
fold    = mod(randperm(max(gid)), cv_folds) + 1;   % assign each bout to a fold
rowfold = fold(gid);

% ── Numeric design for penalised logistic fit ────────────────────────────
% Columns: [intercept | kernel basis | baseline basis | mov dummies | slo dummies].
% Built ONCE; CV subsets rows (keeps dummy columns identical across folds, which
% a per-fold fitglm refit does not guarantee when a fold misses a category).
[~,~,im] = unique(mov);  Dmov = dummyvar(im);  Dmov = Dmov(:,2:end);   % drop reference level
[~,~,is] = unique(slo);  Dslo = dummyvar(is);  Dslo = Dslo(:,2:end);
% Mundlak: append the per-bout mean as an (unpenalised) covariate. The kernel
% (built from within-bout-centred s) then carries the WITHIN-bout effect and this
% column carries the BETWEEN-bout level effect.
if center_per_bout
    Xdes      = [ones(r,1), Xk, Bb, Dmov, Dslo, bmean];
    bmean_col = size(Xdes,2);
else
    Xdes      = [ones(r,1), Xk, Bb, Dmov, Dslo];
    bmean_col = [];
end
ker_cols = 1 + (1:nb_k);                  % kernel-coefficient indices within Xdes
t_rel    = -lag_fr / fps;                 % time relative to escape (s): <0 past, >0 future

% ── Roughness (P-spline) penalty on the kernel coefficients ONLY ─────────
% Raised-cosine bumps are evenly spaced, so a 2nd-difference penalty on the basis
% weights penalises the CURVATURE of β(τ) (Eilers & Marx P-splines). This is what
% damps the +,−,+ anti-correlated solution that overlapping bumps + autocorrelated
% sm produce, without shrinking amplitude per se and without touching the
% baseline/controls. The penalty's 2D null space (constant + linear coeff trends)
% stays unpenalised, so a genuine smooth/flat kernel is unaffected.
switch shuffle_mode
    case 'full'                           % one continuous basis → one penalty
        Dk = diffmat(nb_k, 2);
    case 'split'                          % don't difference across the τ=0 seam
        Dk = blkdiag(diffmat(nb_kernel, 2), diffmat(nb_acausal, 2));
end
Pk = zeros(size(Xdes,2));  Pk(ker_cols,ker_cols) = Dk' * Dk;   % scaled by λ below

% ── Diagnostics: bump overlap + design collinearity ──────────────────────
Rk     = corrcoef(Xk);                    % correlations among kernel design columns
offd   = ~eye(nb_k);
VIF    = diag(inv(Rk));                    % variance-inflation per kernel column
Gb     = Bk' * Bk; dG = sqrt(diag(Gb)); Gcos = Gb ./ (dG*dG');   % cosine sim between bumps
active = sum(Bk > 1e-6, 2);               % how many bumps overlap each lag
fprintf('\n── Kernel-design diagnostics ──\n');
fprintf('design cols: mean|corr|=%.2f  max|corr|=%.2f  cond(corr)=%.0f  max VIF=%.1f\n', ...
    mean(abs(Rk(offd))), max(abs(Rk(offd))), cond(Rk), max(VIF));
fprintf('basis: %d bumps over %.1f s  (mean %.1f bumps active/lag, adjacent-bump cos=%.2f)\n', ...
    nb_k, (max(lag_fr)-min(lag_fr))/fps, mean(active), mean(diag(Gcos,1)));

% ── Select penalty strength λ by grouped CV (1-SE rule) ──────────────────
% The CV curve is typically FLAT across small λ (the oscillations are predictively
% neutral) then drops once the penalty starts flattening the genuine escape peak.
% Plain argmax then sits in the flat region and barely smooths. The 1-SE rule takes
% the LARGEST λ whose CV score is still within 1 SE of the best — the smoothest
% kernel the data do not reject (the standard glmnet default).
if use_penalty
    lam_grid = logspace(-1, 6, 15);            % 0.1 (≈unpenalised) → 1e6 (very smooth)
    cvll_lam = nan(size(lam_grid));  cvse_lam = nan(size(lam_grid));
    for il = 1:numel(lam_grid)
        [cvll_lam(il), llf] = pen_grouped_cv(Xdes, yev, lam_grid(il)*Pk, rowfold, cv_folds);
        cvse_lam(il) = std(llf, 'omitnan') / sqrt(sum(~isnan(llf)));   % SE across folds
    end
    [~, ibest] = max(cvll_lam);
    thr    = cvll_lam(ibest) - cvse_lam(ibest);   % 1 SE below the best score
    i1se   = find(cvll_lam >= thr, 1, 'last');    % largest λ still within 1 SE → most smoothing
    lambda = lam_grid(i1se);
    if i1se == numel(lam_grid)
        fprintf('NOTE: 1-SE λ hit the top grid edge (%.3g) — consider widening lam_grid.\n', lambda);
    end
    fprintf('penalty λ = %.3g  [1-SE rule; argmax was %.3g]  (grid %.2g–%.2g)\n', ...
        lambda, lam_grid(ibest), lam_grid(1), lam_grid(end));
else
    lam_grid = 0; cvll_lam = NaN; cvse_lam = NaN; ibest = 1; i1se = 1; lambda = 0;
    fprintf('penalty DISABLED (use_penalty=false) → unpenalised kernel.\n');
end
P_pen = lambda * Pk;

% ── Fit: unpenalised (fitglm, for the overlay) + penalised (reported) ────
T = table();  T.y = yev;
for j = 1:nb_k,        T.(sprintf('k%d',j)) = Xk(:,j); end
for j = 1:nb_base_use, T.(sprintf('h%d',j)) = Bb(:,j); end
T.mov = categorical(mov);  T.slo = categorical(slo);
kterms = strjoin(arrayfun(@(j) sprintf('k%d',j), 1:nb_k,        'uni', 0), ' + ');
hterms = strjoin(arrayfun(@(j) sprintf('h%d',j), 1:nb_base_use, 'uni', 0), ' + ');
f_full = sprintf('y ~ %s + %s + mov + slo', kterms, hterms);
if center_per_bout, T.bmean = bmean; f_full = [f_full ' + bmean']; end
mdl    = fitglm(T, f_full, 'Distribution','binomial', 'Link','logit');
cn = mdl.CoefficientNames;  kidx = zeros(nb_k,1);
for j = 1:nb_k, kidx(j) = find(strcmp(cn, sprintf('k%d',j))); end
kernel_unpen = Bk * mdl.Coefficients.Estimate(kidx);

[beta_all, V_all] = penalized_logit(Xdes, yev, P_pen);
beta_k = beta_all(ker_cols);
Cov_k  = V_all(ker_cols, ker_cols);
kernel    = Bk * beta_k;                          % per-lag weight (penalised)
kernel_se = sqrt(sum((Bk * Cov_k) .* Bk, 2));     % var = bᵀ Σ b per lag

% ════════════════════════════════════════════════════════════════════════
% Nuisance read-out + ablations
% The baseline(time-in-freeze) and the mov/slo controls are fit but otherwise
% marginalised. Inspect them, and check how much the kernel actually leans on
% them: re-fit dropping controls / baseline / both, holding λ fixed so the ONLY
% thing that changes is the nuisance set. A kernel that barely moves ⇒ the
% nuisances are not confounding it; a big move ⇒ they were carrying real
% variance the kernel would otherwise absorb.
% ════════════════════════════════════════════════════════════════════════
% (a) pull the already-fitted nuisance coefficients out of beta_all / V_all
base_cols = ker_cols(end) + (1:nb_base_use);
nctrl     = size(Dmov,2) + size(Dslo,2);
ctrl_cols = base_cols(end) + (1:nctrl);
beta0     = beta_all(1);
beta_base = beta_all(base_cols);
beta_ctrl = beta_all(ctrl_cols);
se_ctrl   = sqrt(diag(V_all(ctrl_cols, ctrl_cols)));

mov_lv = unique(mov);  slo_lv = unique(slo);      % reference (first) level was dropped
ctrl_lab = [arrayfun(@(v) sprintf('mov=%g', v), mov_lv(2:end), 'uni',0); ...
            arrayfun(@(v) sprintf('slo=%g', v), slo_lv(2:end), 'uni',0)];
fprintf('\n── Nuisance coefficients (log-odds; kernel/controls at reference) ──\n');
fprintf('intercept = %+.3f\n', beta0);
for j = 1:nctrl
    fprintf('  %-10s = %+.3f ± %.3f  (z=%+.1f)\n', ...
        ctrl_lab{j}, beta_ctrl(j), se_ctrl(j), beta_ctrl(j)/max(se_ctrl(j),eps));
end
if nctrl == 0, fprintf('  (no controls: mov/slo each had a single level)\n'); end

% Mundlak between-bout level coefficient (log-odds per 1 SD of bout-mean sm).
% With center_per_bout, the kernel above is WITHIN-bout; this is the BETWEEN-bout
% effect that the simpler (global-z) kernel would otherwise have absorbed.
if center_per_bout && ~isempty(bmean_col)
    beta_bmean = beta_all(bmean_col);
    se_bmean   = sqrt(V_all(bmean_col, bmean_col));
    fprintf('between-bout mean sm (Mundlak) = %+.3f ± %.3f  (z=%+.1f)  [kernel = WITHIN-bout]\n', ...
        beta_bmean, se_bmean, beta_bmean/max(se_bmean,eps));
else
    beta_bmean = NaN; se_bmean = NaN;
    fprintf('per-bout centering OFF → kernel mixes within- and between-bout sm.\n');
end

% (b) baseline hazard vs time-in-freeze (kernel at sm=mean → z=0 → Xk·β=0)
base_logit = beta0 + Bb * beta_base;
[ts, iso]  = sort(tinf);
fh_base = figure('Color','w','Position',[90 90 760 420]);
yyaxis left
plot(ts, 1./(1+exp(-base_logit(iso))), 'LineWidth', 2)
ylabel('baseline hazard  P(unfreeze \mid frozen, sm = mean)')
yyaxis right
histogram(tinf(yev==1), 30, 'FaceAlpha',0.25, 'EdgeColor','none')
ylabel('# observed un-freezes')
xlabel('time in freeze (s)')
title('Fitted baseline hazard vs time-in-freeze')

% (c) ablations: drop controls / baseline / both, same λ, refit the kernel
abl_name = {'full','no controls','no baseline','kernel only'};
abl_des  = { Xdes, [ones(r,1) Xk Bb], [ones(r,1) Xk Dmov Dslo], [ones(r,1) Xk] };
abl_ker  = nan(nLag, numel(abl_des));
for a = 1:numel(abl_des)
    Xa = abl_des{a};
    Pa = zeros(size(Xa,2));  Pa(ker_cols,ker_cols) = Dk' * Dk;   % kernel cols stay at 1+(1:nb_k)
    ba = penalized_logit(Xa, yev, lambda * Pa);
    abl_ker(:,a) = Bk * ba(ker_cols);
end
pk = max(abs(abl_ker(:,1)));
fprintf('\n── Kernel sensitivity to nuisance terms (max |Δβ(τ)| vs full, %% of full peak) ──\n');
for a = 2:numel(abl_des)
    fprintf('  %-12s : %.1f%%\n', abl_name{a}, 100*max(abs(abl_ker(:,a)-abl_ker(:,1)))/max(pk,eps));
end

fh_abl = figure('Color','w','Position',[100 100 760 460]); hold on
cmap_a = lines(numel(abl_des));
for a = 1:numel(abl_des)
    plot(t_rel, abl_ker(:,a), 'Color',cmap_a(a,:), 'LineWidth', 2-0.4*(a>1));
end
yline(0,'k:'); xline(0,'--k','escape')
xlabel('Time rel. to un-freeze (s)'); ylabel('\beta(\tau)')
title('Kernel sensitivity: dropping controls / baseline(t)  (\lambda fixed)')
legend(abl_name, 'Location','best', 'box','off')

% ── Diagnostics figure: why the raw kernel rings, and the penalty's effect ─
fh_diag = figure('Color','w','Position',[60 60 1100 760]);
tiledlayout(2, 2, 'TileSpacing','compact', 'Padding','compact')

nexttile; plot(t_rel, Bk, 'LineWidth', 1); hold on; xline(0,'--k','escape')
xlabel('Time rel. to un-freeze (s)'); ylabel('basis weight')
title(sprintf('Raised-cosine basis: %d bumps, %.1f overlapping/lag', nb_k, mean(active)))

nexttile; imagesc(abs(Rk)); axis square; colorbar; caxis([0 1])
xlabel('kernel basis column'); ylabel('kernel basis column')
title(sprintf('|corr| of design cols  (cond=%.0f, max VIF=%.1f)', cond(Rk), max(VIF)))

nexttile; hold on; grid on
if use_penalty
    errorbar(lam_grid, cvll_lam, cvse_lam, 'o-', 'LineWidth', 1.2, 'CapSize', 3)
    set(gca, 'XScale', 'log')
    yline(cvll_lam(ibest) - cvse_lam(ibest), ':', '1 SE', 'Color',[0.4 0.4 0.4])
    plot(lam_grid(ibest), cvll_lam(ibest), 'ko', 'MarkerSize', 10, 'LineWidth', 1.4)   % argmax
    plot(lambda, cvll_lam(i1se), 'rp', 'MarkerFaceColor','r', 'MarkerSize', 14)         % 1-SE (used)
    legend({'CV \pm SE','1 SE below best','argmax','1-SE (used)'}, 'Location','southwest', 'box','off')
    title('Penalty selection: 1-SE rule (\star used)')
else
    text(0.5, 0.5, 'penalty disabled (unpenalised run)', 'Units','normalized', ...
        'HorizontalAlignment','center')
    title('Penalty selection')
end
xlabel('penalty \lambda'); ylabel('grouped CV logLik/obs')

nexttile; hold on
plot(t_rel, kernel_unpen, 'Color',[0.85 0.4 0.2], 'LineWidth', 1.5)
fill([t_rel; flipud(t_rel)], [kernel+1.96*kernel_se; flipud(kernel-1.96*kernel_se)], ...
    [0.2 0.5 0.8], 'FaceAlpha', 0.2, 'EdgeColor','none')
plot(t_rel, kernel, 'Color',[0.2 0.5 0.8], 'LineWidth', 2)
yline(0,'k:'); xline(0,'--k','escape')
xlabel('Time rel. to un-freeze (s)'); ylabel('\beta(\tau)')
title('Kernel: unpenalised (orange) vs penalised (blue)')
legend({'unpenalised','penalised 95% CI','penalised'}, 'box','off', 'Location','best')

% ── Diagnostic: is the kernel just sm autocorrelation / the event average? ─
% A symmetric β(τ) (causal past mirroring the acausal future) is the fingerprint
% of an autocorrelated regressor: sm is elevated AROUND the un-freeze, and because
% sm's autocorrelation R(τ)=R(−τ) is symmetric, a collinear design returns ≈R(τ)
% rather than a directed/causal kernel. Overlay (shape-normalised) the recovered
% kernel, the sm autocorrelation, and the event-triggered sm excess: if the kernel
% tracks the symmetric R(τ), the fine-lag shape is autocorrelation, not causality.
zero_idx = find(lag_fr == 0, 1);
ev       = yev == 1;
eta_ev   = mean(S(ev,:), 1)';                 % event-triggered average of sm history/future
eta_exc  = eta_ev - mean(S, 1)';              % excess vs the at-risk mean (what marks events)
acf_sm   = mean(S .* S(:,zero_idx), 1)';      % E[s(τ)·s(0)] over all at-risk intervals
acf_sm   = acf_sm / acf_sm(zero_idx);         % normalise so R(0)=1
nrm      = @(v) v / max(abs(v));              % shape-normalise to unit peak for overlay

rho_acf = corr(kernel, acf_sm);
rho_eta = corr(kernel, eta_exc);
fprintf('\n── Kernel-shape diagnostics ──\n');
fprintf('corr(kernel, sm-ACF) = %.2f   corr(kernel, event-triggered excess) = %.2f\n', rho_acf, rho_eta);
fprintf('(high corr with the symmetric sm-ACF ⇒ fine-lag shape is autocorrelation, not directed causality)\n');

fh_acf = figure('Color','w','Position',[80 80 760 460]); hold on
plot(t_rel, nrm(kernel),  'Color',[0.2 0.5 0.8], 'LineWidth', 2)
plot(t_rel, nrm(acf_sm),  'k--', 'LineWidth', 1.6)
plot(t_rel, nrm(eta_exc), 'Color',[0.85 0.4 0.2], 'LineWidth', 1.6)
yline(0,'k:'); xline(0,'--k','escape')
xlabel('Time relative to un-freeze (s)   —   past < 0 < future')
ylabel('shape (normalised to unit peak)')
title(sprintf('Kernel vs sm-autocorrelation / event average  (\\rho_{ACF}=%.2f)', rho_acf))
legend({'recovered kernel','sm autocorr R(\tau)','event-triggered sm excess'}, ...
    'Location','best', 'box','off')

fh_kernel = figure('Color','w','Position',[100 100 720 460]); hold on
xr = [t_rel; flipud(t_rel)];
yr = [kernel + 1.96*kernel_se; flipud(kernel - 1.96*kernel_se)];
fill(xr, yr, [0.2 0.5 0.8], 'FaceAlpha', 0.2, 'EdgeColor','none')
plot(t_rel, kernel, 'Color', [0.2 0.5 0.8], 'LineWidth', 2)
yline(0, 'k:'); xline(0, '--k', 'escape')
% shade the acausal control region
yl = ylim;
patch([0 kernel_future_s kernel_future_s 0], [yl(1) yl(1) yl(2) yl(2)], ...
    [0.9 0.9 0.9], 'FaceAlpha', 0.4, 'EdgeColor','none')
uistack(findobj(gca,'Type','patch'),'bottom')
text(kernel_future_s/2, yl(2)*0.9, 'acausal control', 'HorizontalAlignment','center', 'FontSize',9)
xlabel('Time relative to un-freeze (s)   —   past < 0 < future')
ylabel('Kernel weight  \beta(\tau)   (hazard, per SD of sm)')
title('Social-motion → un-freeze hazard kernel')

% ── Leak time-constant: exponential fit to the causal kernel ─────────────
% Leaky integrator ⇒ causal kernel β(τ) ≈ A·exp(-τ/τ_leak) + c, where τ is the
% lag before escape (s). Weighted by 1/SE² so noisy deep-past lags don't drive
% the fit. CI on τ_leak via parametric bootstrap: resample the basis coeffs
% from the GLM covariance, reconstruct the kernel, refit — no data re-run.
cm     = caus_mask;
tau_s  = lag_fr(cm) / fps;                 % lag before escape (s), ascending from 0
kc     = kernel(cm);
sec    = kernel_se(cm);
w      = 1 ./ max(sec, eps).^2;            % inverse-variance weights

expmod = @(p,t) p(1).*exp(-t./exp(p(2))) + p(3);     % p = [A, log(tau_leak), c]
sse    = @(p) sum(w .* (kc - expmod(p,tau_s)).^2);
p0     = [max(kc)-min(kc), log(0.5), min(kc)];
opt    = optimset('MaxFunEvals',1e4, 'MaxIter',1e4, 'Display','off');
pfit   = fminsearch(sse, p0, opt);
A_hat = pfit(1);  tau_hat = exp(pfit(2));  c_hat = pfit(3);
res    = kc - expmod(pfit, tau_s);
wmean  = sum(w.*kc) / sum(w);
R2_exp = 1 - sum(w.*res.^2) / sum(w.*(kc - wmean).^2);

% model-free cross-check: center-of-mass lag of the positive causal weight
kpos    = max(kc, 0);
tau_com = sum(tau_s .* kpos) / sum(kpos);

% parametric bootstrap CI on tau_leak (resample kernel coeffs ~ N(beta_k, Cov_k))
nboot = 2000;
Csym  = (Cov_k + Cov_k')/2;
[Vc,Dc] = eig(Csym);  Csym = Vc*diag(max(diag(Dc),0))*Vc';   % clip to PSD
Bsamp = mvnrnd(beta_k(:)', Csym, nboot);
Bk_c  = Bk(cm, :);
tau_boot = nan(nboot,1);
for ii = 1:nboot
    kb = Bk_c * Bsamp(ii,:)';
    try
        pb = fminsearch(@(p) sum(w.*(kb - expmod(p,tau_s)).^2), pfit, opt);
        tb = exp(pb(2));
        if tb > 1e-3 && tb < max(tau_s)*5, tau_boot(ii) = tb; end
    catch
    end
end
tau_ci = prctile(tau_boot(~isnan(tau_boot)), [2.5 97.5]);

fprintf('\n── Leak time-constant (exponential fit to causal kernel) ──\n');
fprintf('tau_leak = %.2f s  (95%% CI %.2f–%.2f),  half-life = %.2f s,  weighted R² = %.2f\n', ...
    tau_hat, tau_ci(1), tau_ci(2), tau_hat*log(2), R2_exp);
fprintf('center-of-mass lag (model-free) = %.2f s\n', tau_com);
if R2_exp < 0.5
    fprintf('NOTE: low R² — a single exponential may not describe the kernel well.\n');
end

fh_leak = figure('Color','w','Position',[100 100 560 420]); hold on
fill([tau_s; flipud(tau_s)], [kc+1.96*sec; flipud(kc-1.96*sec)], [0.2 0.5 0.8], ...
    'FaceAlpha',0.15, 'EdgeColor','none')
plot(tau_s, kc, 'Color',[0.2 0.5 0.8], 'LineWidth',2)
plot(tau_s, expmod(pfit, tau_s), 'r--', 'LineWidth',2)
yline(0,'k:')
xlabel('Lag before escape, \tau (s)')
ylabel('Causal kernel  \beta(\tau)')
title(sprintf('Leak fit: \\tau_{leak} = %.2f s (CI %.2f–%.2f), R² = %.2f', ...
    tau_hat, tau_ci(1), tau_ci(2), R2_exp))
legend({'real 95% CI','real kernel','A e^{-\tau/\tau_{leak}}+c'}, 'box','off')

% ── Stage 2: integration vs extrema, by AIC + grouped CV ─────────────────
zc = @(x) (x - mean(x)) ./ std(x);
T.f_mean  = zc(fmean);
T.f_max   = zc(fmax);
T.f_slope = zc(fslp);
ctrl = sprintf('%s + mov + slo', hterms);   % shared baseline + controls

specs = struct( ...
  'name', {'baseline','mean (integration)','max (extrema)','slope (transient)', ...
           'mean+max','mean+slope','mean+max+slope'}, ...
  'rhs',  {ctrl, ...
           ['f_mean + '  ctrl], ...
           ['f_max + '   ctrl], ...
           ['f_slope + ' ctrl], ...
           ['f_mean + f_max + ' ctrl], ...
           ['f_mean + f_slope + ' ctrl], ...
           ['f_mean + f_max + f_slope + ' ctrl]});

% (CV folds `rowfold` already assigned above, before the kernel fit, so the
%  λ-selection, this Stage-2 comparison, and the shuffle null all share folds.)

fprintf('\n%-22s %12s %14s\n','model','AIC','CV logLik/obs');
aic = nan(numel(specs),1); cvll = nan(numel(specs),1);
for m = 1:numel(specs)
    mm      = fitglm(T, ['y ~ ' specs(m).rhs], 'Distribution','binomial');
    aic(m)  = mm.ModelCriterion.AIC;
    cvll(m) = grouped_cv_loglik(T, ['y ~ ' specs(m).rhs], rowfold, cv_folds);
    fprintf('%-22s %12.1f %14.4f\n', specs(m).name, aic(m), cvll(m));
end
[~, best] = max(cvll);
fprintf('\nBest by CV log-likelihood: %s\n', specs(best).name);
fprintf(['Read-out (compare the COMBINED rows to mean-alone for UNIQUE info):\n' ...
         '  mean ≫ baseline                 ⇒ integration/accumulation\n' ...
         '  mean+max  beats mean (ΔAIC>~4)   ⇒ extra peak/extrema sensitivity\n' ...
         '  mean+slope beats mean            ⇒ genuine temporal-derivative (timing) sensitivity\n' ...
         '  max/slope alone but no gain over mean ⇒ redundant with the level.\n']);

% ── Control: time-shuffle of social motion ───────────────────────────────
% 'full'  → permute ALL lags within each row (preserves the window value
%           multiset, destroys all order; flat leakage-free null, no seam).
% 'split' → permute the causal and acausal halves independently (preserves
%           each half's mean; assesses causal vs acausal ordering separately).
% real ≫ shuffled ⇒ temporal ORDER carries info beyond the value distribution.
caus_idx  = find(caus_mask);   nca = numel(caus_idx);
acaus_idx = find(~caus_mask);  naa = numel(acaus_idx);

ll_real = pen_grouped_cv(Xdes, yev, P_pen, rowfold, cv_folds);

ll_null    = nan(n_shuffle, 1);
kernels_sh = nan(nLag, n_shuffle);           % shuffled kernels → flat null band for the SHAPE
nr      = size(S, 1);
rows_all = repmat((1:nr)', 1, nLag);
Ac = S(:, caus_idx);   rows_c = repmat((1:nr)', 1, nca);
Aa = S(:, acaus_idx);  rows_a = repmat((1:nr)', 1, naa);
Xdes_sh = Xdes;                              % template; only the kernel cols get re-shuffled
tic_shuf = tic;
for sN = 1:n_shuffle
    if mod(sN, max(1, round(n_shuffle/20))) == 0
        fprintf('  shuffle %d/%d  (%.1f min elapsed)\n', sN, n_shuffle, toc(tic_shuf)/60);
    end
    switch shuffle_mode
        case 'full'                              % permute all lags within each row
            [~, ord] = sort(rand(nr, nLag), 2);
            Ssh = S(sub2ind([nr nLag], rows_all, ord));
        case 'split'                             % permute each half independently
            [~, ordc] = sort(rand(nr, nca), 2);
            [~, orda] = sort(rand(nr, naa), 2);
            Ssh = S;
            Ssh(:, caus_idx)  = Ac(sub2ind([nr nca], rows_c, ordc));
            Ssh(:, acaus_idx) = Aa(sub2ind([nr naa], rows_a, orda));
    end
    Xdes_sh(:, ker_cols) = Ssh * Bk;         % swap in shuffled kernel projections
    ll_null(sN)      = pen_grouped_cv(Xdes_sh, yev, P_pen, rowfold, cv_folds);
    b_sh             = penalized_logit(Xdes_sh, yev, P_pen);
    kernels_sh(:,sN) = Bk * b_sh(ker_cols);  % full-data penalised refit (same λ)
end

pval = (1 + sum(ll_null >= ll_real)) / (1 + n_shuffle);
fprintf('\n── Time-shuffle control (mode = %s) ──\n', shuffle_mode);
fprintf('Real kernel     CV logLik/obs : %.5f\n', ll_real);
fprintf('Shuffled (n=%d)               : %.5f ± %.5f  [min %.5f, max %.5f]\n', ...
    n_shuffle, mean(ll_null), std(ll_null), min(ll_null), max(ll_null));
fprintf('p(shuffled ≥ real)            : %.3f\n', pval);
% sanity: the mean preserved by the shuffle should be perfectly correlated (≈1)
switch shuffle_mode
    case 'full'
        fprintf('sanity  corr(window-mean, shuffled) = %.3f (expect 1)\n', ...
            corr(mean(S,2), mean(Ssh,2)));
    case 'split'
        fprintf('sanity  corr(causal-mean, shuffled) = %.3f, corr(acausal-mean, shuffled) = %.3f (expect 1)\n', ...
            corr(mean(Ac,2), mean(Ssh(:,caus_idx),2)), corr(mean(Aa,2), mean(Ssh(:,acaus_idx),2)));
end
fprintf(['Interpretation: real ≫ shuffled ⇒ temporal ORDER matters (kernel\n' ...
         'shape is real). real ≈ shuffled ⇒ only the value distribution\n' ...
         '(level/peaks) matters, not when within the window it occurs.\n']);

% shuffled-kernel null band: parametric mean ± 1.96 SD across shuffles.
% Smoother than empirical percentiles at fixed n_shuffle, and the mean line
% shows the (flat) null center directly. Real kernel exiting the band ⇒ timing.
sh_mu = mean(kernels_sh, 2);
sh_sd = std(kernels_sh, 0, 2);
sh_lo = sh_mu - 1.96*sh_sd;
sh_hi = sh_mu + 1.96*sh_sd;

fh_ctrl = figure('Color','w','Position',[100 100 1040 420]);
tiledlayout(1, 2, 'TileSpacing','compact', 'Padding','compact')

nexttile; hold on   % real kernel vs shuffled-kernel null band
fill([t_rel; flipud(t_rel)], [sh_hi; flipud(sh_lo)], [0.6 0.6 0.6], ...
    'FaceAlpha', 0.35, 'EdgeColor','none')
plot(t_rel, sh_mu, 'Color', [0.4 0.4 0.4], 'LineWidth', 1)   % flat null center
fill([t_rel; flipud(t_rel)], [kernel + 1.96*kernel_se; flipud(kernel - 1.96*kernel_se)], ...
    [0.2 0.5 0.8], 'FaceAlpha', 0.2, 'EdgeColor','none')
plot(t_rel, kernel, 'Color', [0.2 0.5 0.8], 'LineWidth', 2)
yline(0, 'k:'); xline(0, '--k', 'escape')
xlabel('Time relative to un-freeze (s)   —   past < 0 < future')
ylabel('Kernel weight  \beta(\tau)')
title('Real kernel vs time-shuffled null (mean \pm 1.96 SD)')
legend({'shuffled null (95%)','null mean','real 95% CI','real kernel'}, 'Location','best', 'box','off')

nexttile; hold on   % CV log-likelihood: real vs shuffled null
histogram(ll_null, max(5, round(n_shuffle/3)), 'FaceColor', [0.7 0.7 0.7], 'EdgeColor','none')
yl = ylim;
plot([ll_real ll_real], yl, 'r', 'LineWidth', 2)
text(ll_real, yl(2)*0.95, sprintf('  real  (p=%.3f)', pval), 'Color','r', 'FontWeight','bold')
xlabel('CV log-likelihood / obs')
ylabel('Count (shuffled nulls)')
title('Predictive power: real vs shuffle')
legend({'shuffled null','real'}, 'Location','northwest', 'box','off')

% ── Persist results + figures (overnight runs) ───────────────────────────
if save_out
    stamp  = datestr(now, 'yyyymmdd_HHMMSS'); %#ok<TNOW1,DATST>
    outdir = paths.dataset;
    try
        params = struct('kernel_past_s',kernel_past_s, 'kernel_future_s',kernel_future_s, ...
            'dt_frames',dt_frames, 'nb_kernel',nb_kernel, 'nb_acausal',nb_acausal, ...
            'bump_width',bump_width, 'nb_baseline',nb_baseline, 'cv_folds',cv_folds, ...
            'n_shuffle',n_shuffle, 'contact_threshold',contact_threshold, 'trunc_point',trunc_point, ...
            'shuffle_mode',shuffle_mode, 'use_penalty',use_penalty, 'lambda',lambda, ...
            'anchor_zero',anchor_zero, 'center_per_bout',center_per_bout);
        save(fullfile(outdir, ['hazard_kernel_' stamp '.mat']), ...
            'kernel','kernel_se','kernel_unpen','t_rel','lag_fr','kernels_sh','sh_mu','sh_sd', ...
            'll_real','ll_null','pval','aic','cvll','specs','params', ...
            'tau_hat','tau_ci','R2_exp','tau_com','A_hat','c_hat','tau_boot', ...
            'beta_k','Cov_k','Bk','caus_mask','fps', 'use_penalty', ...
            'lambda','lam_grid','cvll_lam','cvse_lam','Rk','VIF','active', ...
            'acf_sm','eta_ev','eta_exc','rho_acf','rho_eta', ...
            'beta0','beta_base','beta_ctrl','se_ctrl','ctrl_lab','base_logit', ...
            'abl_ker','abl_name','center_per_bout','beta_bmean','se_bmean', '-v7.3');
        exportgraphics(fh_kernel, fullfile(outdir, ['hazard_kernel_'  stamp '.png']), 'Resolution',200);
        exportgraphics(fh_ctrl,   fullfile(outdir, ['hazard_control_' stamp '.png']), 'Resolution',200);
        exportgraphics(fh_leak,   fullfile(outdir, ['hazard_leakfit_' stamp '.png']), 'Resolution',200);
        exportgraphics(fh_diag,   fullfile(outdir, ['hazard_diag_'    stamp '.png']), 'Resolution',200);
        exportgraphics(fh_acf,    fullfile(outdir, ['hazard_acf_'     stamp '.png']), 'Resolution',200);
        exportgraphics(fh_base,   fullfile(outdir, ['hazard_baseline_' stamp '.png']), 'Resolution',200);
        exportgraphics(fh_abl,    fullfile(outdir, ['hazard_ablation_' stamp '.png']), 'Resolution',200);
        fprintf('\nSaved results + figures to %s (stamp %s)\n', outdir, stamp);
    catch ME
        warning('save_out failed: %s', ME.message);
    end
end

% ════════════════════════════════════════════════════════════════════════
% Local functions
% ════════════════════════════════════════════════════════════════════════
function B = raised_cosine_basis(x, nb, width_mult, anchor)
% Linearly-spaced raised-cosine bumps tiling the range of x.
% width_mult sets each bump's half-width in units of the center spacing:
%   1   → adjacent bumps only (minimal overlap; sums to ≈1 but per-lag variance
%         scallops with the spacing → bumpy reconstruction band)
%   2-3 → 3-4 bumps overlap everywhere → near-uniform leverage → smooth band.
% anchor (optional): if given, the whole evenly-spaced grid is shifted by < half a
%   spacing so one center lands exactly on `anchor` (e.g. 0 → a knot at τ=0). Spacing
%   stays uniform (the 2nd-difference penalty is unaffected); edge bumps shift slightly.
    if nargin < 3 || isempty(width_mult), width_mult = 1; end
    if nargin < 4, anchor = []; end
    x = x(:);
    lo = min(x); hi = max(x);
    c  = linspace(lo, hi, nb);
    sp = (hi - lo) / (nb - 1);
    if ~isempty(anchor)
        [~, jn] = min(abs(c - anchor));   % center nearest the anchor
        c = c + (anchor - c(jn));         % shift grid so that center sits exactly on anchor
    end
    w  = width_mult * sp;          % bump half-width (support ±w)
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
% Penalised logistic regression by IRLS / Newton:
%   maximise  loglik(y | Xb)  −  0.5 * bᵀ P b.
% P is the (p×p) quadratic penalty (zeros for unpenalised coefficients). V is the
% penalised posterior covariance inv(XᵀWX + P) — the P-spline Bayesian covariance,
% used for the kernel CI (it correctly narrows where the penalty pins the shape).
    [~, p] = size(X);
    b   = zeros(p, 1);
    rdg = 1e-8 * eye(p);                 % tiny jitter for numerical stability
    for it = 1:100 %#ok<NASGU>
        eta = min(max(X*b, -30), 30);
        mu  = 1 ./ (1 + exp(-eta));
        w   = max(mu .* (1 - mu), 1e-9);
        H   = X' * (X .* w) + P + rdg;   % penalised Hessian = XᵀWX + P
        g   = X' * (y - mu) - P * b;     % penalised gradient
        step = H \ g;
        b   = b + step;
        if max(abs(step)) < 1e-8, break; end
    end
    V = inv(H);
end

function [ll, ll_folds] = pen_grouped_cv(X, y, P, rowfold, K)
% Held-out Bernoulli log-likelihood/obs of the penalised model, folds grouped by
% bout (rowfold). Same fixed design X across folds; only rows are subset.
% Returns the across-fold mean `ll` and the per-fold means `ll_folds` (the latter
% feeds the 1-SE rule's standard error).
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

function ll = grouped_cv_loglik(T, formula, rowfold, K)
% Mean held-out Bernoulli log-likelihood, folds grouped by bout.
    n = height(T); acc = 0; cnt = 0;
    for k = 1:K
        te = rowfold == k; tr = ~te;
        if ~any(te) || ~any(tr), continue; end
        mm = fitglm(T(tr,:), formula, 'Distribution','binomial');
        p  = predict(mm, T(te,:));
        p  = min(max(p, 1e-9), 1-1e-9);
        y  = T.y(te);
        acc = acc + sum(y.*log(p) + (1-y).*log(1-p));
        cnt = cnt + sum(te);
    end
    ll = acc / cnt;
end
