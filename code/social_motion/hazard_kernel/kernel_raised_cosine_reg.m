clearvars; clc;
% ════════════════════════════════════════════════════════════════════════
% Raised-cosine hazard kernel on the REAL data: WITH vs WITHOUT smoothness.
%
% Fits the social-motion → un-freeze kernel β(τ) in the raised-cosine basis two
% ways, on the same design used by hazard_kernel_sparse.m:
%   • UNPENALISED   — plain logistic regression on the bump projections.
%   • SMOOTHNESS    — P-spline (2nd-difference) roughness penalty on the bump
%                     coefficients, penalty strength λ chosen by grouped CV
%                     (1-SE rule). This is the paper's Table-1 "smoothness" prior
%                     and the method in hazard_kernel_sm.m.
% The overlay shows what the smoothness penalty buys: the unpenalised kernel
% typically RINGS (overlapping bumps + autocorrelated sm → +,−,+ oscillation);
% the penalised one damps that curvature.
% ════════════════════════════════════════════════════════════════════════

fps             = 60;
kernel_past_s   = 9.0;
kernel_future_s = 3.0;
dt_frames       = 3;
nb_kernel       = 27;
nb_acausal      = 9;
bump_width      = 1.25;
nb_baseline     = 6;
cv_folds        = 10;
trunc_point     = 30;   entry_fr = trunc_point;
mask_preonset   = false; % true → zero out lags reaching before freeze onset
grid_anchor     = 'entry';% 'entry' = forward from onset | 'offset' = backward from the freeze break
contact_threshold = 70;
save_out        = true;
rng(10);

% ── Load + filter (same pipeline as hazard_kernel_sparse.m) ───────────────
threshold_imm = 2; threshold_mob = 2; threshold_pc = 4;
id_code = sprintf('imm%d_mob%d_pc%d', threshold_imm, threshold_mob, threshold_pc);
paths        = path_generator('folder','social_motion/hazard_kernel','bouts_id',id_code,'imfirst',false);
bouts        = importdata(fullfile(paths.dataset,  'bouts.mat'));
motion_cache = importdata(fullfile(paths.cache_path,'motion_cache.mat'));
thresholds = define_thresholds('le_window', struct('le_window_sl',[5 55],'le_window_fl',[5 55]));
bouts = bouts_formatting(bouts, thresholds);
bl = data_parser_new(bouts, 'type','immobility','period','loom','window','le','nloom',2:20,'min_dur',trunc_point);
bl = impose_contact_threshold(bl, 'threshold', contact_threshold, 'type','onlyfreeze');
bl.dur_frames = bl.ending_time;
bl.bout_id    = (1:height(bl))';
bl = bl(bl.dur_frames > trunc_point, :);
fprintf('Bouts: %d total, %d censored (%.1f%%)\n', height(bl), sum(bl.is_censored), 100*mean(bl.is_censored));

% ── Lag grid ──────────────────────────────────────────────────────────────
back_fr = round(kernel_past_s*fps);  fwd_fr = round(kernel_future_s*fps);
lag_fr  = (-fwd_fr:back_fr)';  nLag = numel(lag_fr);
t_rel   = -lag_fr / fps;

% ── Global z-scoring of sm ────────────────────────────────────────────────
flies = unique(bl.fly);  nn=0; s1=0; s2=0;
for i = 1:numel(flies)
    v = motion_cache(flies(i)); v = v(:); v = v(~isnan(v));
    nn = nn + numel(v); s1 = s1 + sum(v); s2 = s2 + sum(v.^2);
end
sm_mu = s1/nn;  sm_sd = sqrt(s2/nn - sm_mu^2);

% ── Build long (person-period) design ─────────────────────────────────────
max_rows = sum(ceil(bl.dur_frames / dt_frames)) + height(bl);
S=nan(max_rows,nLag); tinf=nan(max_rows,1); yev=zeros(max_rows,1);
gid=nan(max_rows,1); mov=nan(max_rows,1); slo=nan(max_rows,1);
r = 0;
for b = 1:height(bl)
    sm = motion_cache(bl.fly(b)); sm = sm(:);  L = numel(sm);
    on = bl.onsets(b);  dur = bl.dur_frames(b);
    if dur < entry_fr, continue; end
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
        f = grid(gi);  t_abs = on + f - 1;
        if (t_abs - back_fr) < 1 || (t_abs + fwd_fr) > L, continue; end
        s = (sm(t_abs - lag_fr) - sm_mu) / sm_sd;
        if any(isnan(s)), continue; end
        if mask_preonset, s(lag_fr > f-1) = 0; end   % drop pre-onset drive
        r = r + 1;
        S(r,:) = s'; tinf(r) = f/fps; gid(r) = bl.bout_id(b);
        mov(r) = bl.moving_flies(b); slo(r) = bl.sloom(b);
        if gi == numel(grid) && ~bl.is_censored(b), yev(r) = 1; end
    end
end
S=S(1:r,:); tinf=tinf(1:r); yev=yev(1:r); gid=gid(1:r); mov=mov(1:r); slo=slo(1:r);
fprintf('Design: %d at-risk intervals, %d un-freeze events (hazard %.3f)\n', r, sum(yev), mean(yev));

% ── Nuisance block: baseline(time-in-freeze) + mov/slo dummies ────────────
ut = tiedrank(tinf) / numel(tinf);
Bb = raised_cosine_basis(ut, nb_baseline);  Bb = Bb(:, 2:end);
[~,~,im]  = unique(mov);  Dmov = dummyvar(im);  Dmov = Dmov(:,2:end);
[~,~,iss] = unique(slo);  Dslo = dummyvar(iss); Dslo = Dslo(:,2:end);
Nuis = [Bb, Dmov, Dslo];

fold = mod(randperm(max(gid)), cv_folds) + 1;  rowfold = fold(gid);

% ── Raised-cosine kernel basis + design ───────────────────────────────────
Bk    = raised_cosine_basis(lag_fr, nb_kernel + nb_acausal, bump_width);
nb_k  = size(Bk,2);
X     = [ones(r,1), S*Bk, Nuis];
kcols = 1 + (1:nb_k);
Dk    = diffmat(nb_k, 2);
Pk    = zeros(size(X,2));  Pk(kcols,kcols) = Dk'*Dk;   % scaled by λ below

% ── (0) nuisance-only baseline (reference) ────────────────────────────────
Xnull = [ones(r,1), Nuis];  Ptiny = 1e-6*eye(size(Xnull,2));
[bN, VN] = penalized_logit(Xnull, yev, Ptiny);
[devN, dfN, aicN] = model_metrics_l2(Xnull, yev, bN, VN);
cv_devN = -2 * pen_grouped_cv(Xnull, yev, Ptiny, rowfold, cv_folds);

% ── (1) UNPENALISED raised-cosine kernel ──────────────────────────────────
[bU, VU] = penalized_logit(X, yev, zeros(size(X,2)));
kernel_unpen    = Bk * bU(kcols);
kernel_unpen_se = sqrt(sum((Bk * VU(kcols,kcols)) .* Bk, 2));
[devU, dfU, aicU] = model_metrics_l2(X, yev, bU, VU);
cv_devU = -2 * pen_grouped_cv(X, yev, zeros(size(X,2)), rowfold, cv_folds);

% ── (2) SMOOTHNESS-PENALISED raised-cosine kernel (λ by CV, 1-SE) ─────────
lam_grid = logspace(-3, 8, 30);
[lambda, cvllP, hit_edge] = select_lambda_l2(X, yev, Pk, lam_grid, rowfold, cv_folds);
[bP, VP] = penalized_logit(X, yev, lambda*Pk);
kernel_pen    = Bk * bP(kcols);
kernel_pen_se = sqrt(sum((Bk * VP(kcols,kcols)) .* Bk, 2));
[devP, dfP, aicP] = model_metrics_l2(X, yev, bP, VP);
cv_devP = -2 * cvllP;

% ── Report ────────────────────────────────────────────────────────────────
fprintf('\n── Raised-cosine kernel: unpenalised vs smoothness (lower AIC / CVdev = better) ──\n');
fprintf('%-24s %8s %6s %10s %12s\n','model','dev','df','AIC','CVdev/obs');
fprintf('%-24s %8.1f %6.1f %10.1f %12.4f\n','(0) nuisance-only', devN, dfN, aicN, cv_devN);
fprintf('%-24s %8.1f %6.1f %10.1f %12.4f\n','(1) unpenalised',   devU, dfU, aicU, cv_devU);
fprintf('%-24s %8.1f %6.1f %10.1f %12.4f\n','(2) smoothness',    devP, dfP, aicP, cv_devP);
fprintf('smoothness λ = %.3g  (grid %.2g–%.2g)%s\n', lambda, lam_grid(1), lam_grid(end), ...
    repmat('  [NOTE: λ hit top grid edge — widen lam_grid]', 1, hit_edge));

% ── Figure: unpenalised vs penalised kernel ───────────────────────────────
fh = figure('Color','w','Position',[100 100 780 500]); hold on
fill([t_rel; flipud(t_rel)], [kernel_pen+1.96*kernel_pen_se; flipud(kernel_pen-1.96*kernel_pen_se)], ...
   [0.2 0.5 0.8], 'FaceAlpha',0.18, 'EdgeColor','none');
pU = plot(t_rel, kernel_unpen, 'Color',[0.85 0.4 0.2], 'LineWidth',2);
pP = plot(t_rel, kernel_pen,   'Color',[0.2 0.5 0.8], 'LineWidth',2.2);
yline(0,'k:'); xline(0,'--k')
yl = ylim;
patch([0 kernel_future_s kernel_future_s 0], [yl(1) yl(1) yl(2) yl(2)], [0.9 0.9 0.9], 'FaceAlpha',0.35, 'EdgeColor','none')
uistack(findobj(gca,'Type','patch'),'bottom')
text(kernel_future_s/2, yl(2)*0.9, 'acausal control', 'HorizontalAlignment','center', 'FontSize',9)
xlabel('time relative to un-freeze (s)')
ylabel('kernel weight  \beta(\tau)')
title('Raised-cosine hazard kernel: unpenalised vs smoothness-regularised')
legend([pU pP], {sprintf('unpenalised (df=%.0f)',dfU), sprintf('smoothness, 95%% CI (df=%.1f)',dfP)}, ...
   'Location','best', 'box','off')
apply_generic(gca)

% ── Persist ───────────────────────────────────────────────────────────────
if save_out
    stamp  = datestr(now,'yyyymmdd_HHMMSS'); %#ok<TNOW1,DATST>
    outdir = paths.dataset;
    try
        save(fullfile(outdir, ['kernel_rc_reg_' stamp '.mat']), ...
            'kernel_unpen','kernel_unpen_se','kernel_pen','kernel_pen_se','t_rel','lag_fr', ...
            'lambda','lam_grid','devU','dfU','aicU','cv_devU','devP','dfP','aicP','cv_devP', ...
            'devN','dfN','aicN','cv_devN','nb_kernel','nb_acausal','bump_width', '-v7.3');
        exportgraphics(fh, fullfile(outdir, ['kernel_rc_reg_' stamp '.png']), 'Resolution',200);
        fprintf('\nSaved to %s (stamp %s)\n', outdir, stamp);
    catch ME
        warning('save_out failed: %s', ME.message);
    end
end

% ════════════════════════════════════════════════════════════════════════
% Local functions
% ════════════════════════════════════════════════════════════════════════
function B = raised_cosine_basis(x, nb, width_mult, anchor)
    if nargin < 3 || isempty(width_mult), width_mult = 1; end
    if nargin < 4, anchor = []; end
    x = x(:);  lo = min(x); hi = max(x);
    c  = linspace(lo, hi, nb);  sp = (hi - lo) / (nb - 1);
    if ~isempty(anchor), [~,jn] = min(abs(c - anchor));  c = c + (anchor - c(jn)); end
    w  = width_mult * sp;  B = zeros(numel(x), nb);
    for j = 1:nb
        d = (x - c(j)) / w;  in = abs(d) <= 1;
        B(in, j) = 0.5 * (1 + cos(pi * d(in)));
    end
end

function D = diffmat(n, ord)
    D = eye(n);  for k = 1:ord, D = diff(D); end
end

function [b, V] = penalized_logit(X, y, P)
    [~, p] = size(X);  b = zeros(p,1);  rdg = 1e-8*eye(p);
    for it = 1:100 %#ok<NASGU>
        eta = min(max(X*b,-30),30);  mu = 1./(1+exp(-eta));
        w = max(mu.*(1-mu), 1e-9);
        H = X'*(X.*w) + P + rdg;  g = X'*(y-mu) - P*b;
        step = H\g;  b = b + step;
        if max(abs(step)) < 1e-8, break; end
    end
    V = inv(H);
end

function [dev, df, aic] = model_metrics_l2(X, y, b, V)
    eta = min(max(X*b,-30),30);  mu = min(max(1./(1+exp(-eta)),1e-9),1-1e-9);
    dev = -2*sum(y.*log(mu) + (1-y).*log(1-mu));
    w   = max(mu.*(1-mu), 1e-9);
    df  = trace(V * (X'*(X.*w)));
    aic = dev + 2*df;
end

function [ll, ll_folds] = pen_grouped_cv(X, y, P, rowfold, K)
    ll_folds = nan(K,1);
    for k = 1:K
        te = rowfold==k; tr = ~te;
        if ~any(te) || ~any(tr), continue; end
        b = penalized_logit(X(tr,:), y(tr), P);
        eta = min(max(X(te,:)*b,-30),30);  mu = min(max(1./(1+exp(-eta)),1e-9),1-1e-9);
        yy = y(te);  ll_folds(k) = mean(yy.*log(mu) + (1-yy).*log(1-mu));
    end
    ll = mean(ll_folds, 'omitnan');
end

function [lambda, cvll, hit_edge] = select_lambda_l2(X, y, Pk, lam_grid, rowfold, K)
    cvll_lam = nan(size(lam_grid));  cvse_lam = nan(size(lam_grid));
    for il = 1:numel(lam_grid)
        [cvll_lam(il), llf] = pen_grouped_cv(X, y, lam_grid(il)*Pk, rowfold, K);
        cvse_lam(il) = std(llf,'omitnan') / sqrt(sum(~isnan(llf)));
    end
    [~, ibest] = max(cvll_lam);
    thr  = cvll_lam(ibest) - cvse_lam(ibest);
    i1se = find(cvll_lam >= thr, 1, 'last');      % largest λ within 1 SE (most smoothing)
    lambda = lam_grid(i1se);  cvll = cvll_lam(i1se);
    hit_edge = (i1se == numel(lam_grid));
end
