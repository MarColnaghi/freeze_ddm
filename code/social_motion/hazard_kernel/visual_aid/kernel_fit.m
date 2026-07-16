clearvars; clc;
% ════════════════════════════════════════════════════════════════════════
% kernel_fit — VISUAL AID  (closes the loop after basis_projection.m)
%
% How we get from the design matrix to the kernel β(τ). Three steps, one figure:
%
%   [ Xₖ | y ]  --fit-->  b  --β=B·b-->  β(τ)
%
%   • Xₖ = S·B is the projected design (as in basis_projection.m); y is the 0/1
%     outcome (1 = the freeze broke this interval). LEFT panel shows a few
%     example at-risk intervals + their outcome.
%   • Logistic regression (IRLS) finds the coefficient on each basis column. We
%     fit it TWO ways: penalised (a smoothness penalty whose strength is chosen by
%     grouped cross-validation — the real hazard_kernel_sm/-sparse method) and
%     unpenalised (ML in the bump basis). MIDDLE panel: the fitted coefficients b,
%     grouped bars comparing the two fits.
%   • The kernel is the reverse of the projection: β(τ) = B·b. RIGHT panel: both
%     kernels (penalised vs unpenalised) with their 95% CIs.
%
% The fit uses ALL at-risk intervals; the LEFT panel just displays a slice.
% ════════════════════════════════════════════════════════════════════════

% ── parameters (mirror hazard_kernel_sparse.m) ───────────────────────────────
fps             = 60;
kernel_past_s   = 4.0;
kernel_future_s = 2;
dt_frames       = 6;
trunc_point     = 30;   entry_fr = trunc_point;
contact_threshold = 70;
grid_anchor     = 'offset';
mask_preonset   = true;

nb_basis        = 12;   % raised-cosine bumps over the lag axis (the kernel basis B)
bump_width      = 1.25;
nb_baseline     = 6;    % baseline-hazard bumps over time-in-freeze (nuisance)
cv_folds        = 5;    % grouped CV folds to pick the smoothness penalty
include_baseline = false;  % false → drop the time-in-freeze baseline hazard block
include_nuisance = false;  % false → drop the moving_flies + loom-speed dummies
save_out        = true;
rng(124);

% ── load + filter (identical to the real pipeline) ───────────────────────────
threshold_imm = 2; threshold_mob = 2; threshold_pc = 4;
id_code = sprintf('imm%d_mob%d_pc%d', threshold_imm, threshold_mob, threshold_pc);
paths        = path_generator('folder','social_motion/hazard_kernel/visual_aid','bouts_id',id_code,'imfirst',false);
bouts        = importdata(fullfile(paths.dataset,  'bouts.mat'));
motion_cache = importdata(fullfile(paths.cache_path,'motion_cache.mat'));
thresholds = define_thresholds('le_window', struct('le_window_sl',[5 55],'le_window_fl',[5 55]));
bouts = bouts_formatting(bouts, thresholds);

bl = data_parser_new(bouts, 'type','immobility','period','loom','window','le','nloom',2:24,'min_dur',trunc_point);
bl = impose_contact_threshold(bl, 'threshold', contact_threshold, 'type','onlyfreeze');
bl.dur_frames = bl.ending_time;
bl.bout_id    = (1:height(bl))';
bl = bl(bl.dur_frames > trunc_point, :);
fprintf('Bouts: %d (%.1f%% censored)\n', height(bl), 100*mean(bl.is_censored));

% ── global z-scoring of sm ───────────────────────────────────────────────────
flies = unique(bl.fly); nn=0; s1=0; s2=0;
for i = 1:numel(flies)
    v = motion_cache(flies(i)); v = v(:); v = v(~isnan(v));
    nn=nn+numel(v); s1=s1+sum(v); s2=s2+sum(v.^2);
end
sm_mu = s1/nn; sm_sd = sqrt(s2/nn - sm_mu^2);

% ── lag grid ─────────────────────────────────────────────────────────────────
back_fr = round(kernel_past_s*fps);  fwd_fr = round(kernel_future_s*fps);
lag_fr  = (-fwd_fr:back_fr)';  nLag = numel(lag_fr);  lag_s = lag_fr/fps;
t_rel   = -lag_s;

% ── build the FULL person-period design S, outcome yev, nuisances ────────────
max_rows = sum(ceil(bl.dur_frames / dt_frames)) + height(bl);
S=nan(max_rows,nLag); yev=zeros(max_rows,1); tinf=nan(max_rows,1);
gid=nan(max_rows,1); mov=nan(max_rows,1); slo=nan(max_rows,1);
r = 0;
for b = 1:height(bl)
    sm  = motion_cache(bl.fly(b)); sm = sm(:); L = numel(sm);
    on  = bl.onsets(b); dur = bl.dur_frames(b);
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
        f = grid(gi); t_abs = on + f - 1;
        if (t_abs - back_fr) < 1 || (t_abs + fwd_fr) > L, continue; end
        s = (sm(t_abs - lag_fr) - sm_mu) / sm_sd;
        if any(isnan(s)), continue; end
        if mask_preonset, s(lag_fr > f-1) = 0; end
        r = r + 1;
        S(r,:) = s'; tinf(r) = f/fps; gid(r) = bl.bout_id(b);
        mov(r) = bl.moving_flies(b); slo(r) = bl.sloom(b);
        if gi == numel(grid) && ~bl.is_censored(b), yev(r) = 1; end
    end
end
S=S(1:r,:); yev=yev(1:r); tinf=tinf(1:r); gid=gid(1:r); mov=mov(1:r); slo=slo(1:r);
fprintf('Design: %d at-risk intervals, %d un-freeze events (hazard %.3f)\n', r, sum(yev), mean(yev));

% ── basis, projection, nuisances ─────────────────────────────────────────────
B  = raised_cosine_basis(lag_fr, nb_basis, bump_width);   % nLag × K
K  = size(B,2);
Xk = S * B;

Nuis = [];
if include_baseline
    ut = tiedrank(tinf) / numel(tinf);
    Bb = raised_cosine_basis(ut, nb_baseline); Bb = Bb(:,2:end);
    Nuis = [Nuis, Bb];
end
if include_nuisance
    [~,~,im] = unique(mov);  Dmov = dummyvar(im);  Dmov = Dmov(:,2:end);
    [~,~,iv] = unique(slo);  Dslo = dummyvar(iv);  Dslo = Dslo(:,2:end);
    Nuis = [Nuis, Dmov, Dslo];
end
fprintf('Design cols: intercept + %d kernel + %d nuisance (baseline=%d, mov/slo=%d)\n', ...
    K, size(Nuis,2), include_baseline, include_nuisance);

% ── fit BOTH: penalised (CV-tuned smoothness) and unpenalised ML ─────────────
X   = [ones(r,1), Xk, Nuis];
ker = 1 + (1:K);                                  % kernel coefficient columns
Dk  = diffmat(K, 2);                              % 2nd-difference (curvature) operator
Pk  = zeros(size(X,2)); Pk(ker,ker) = Dk' * Dk;   % penalty only on the kernel block

% (a) penalised: pick the smoothness strength by grouped CV, then fit on all data
fold    = mod(randperm(max(gid)), cv_folds) + 1;
rowfold = fold(gid);
lam_L2  = logspace(-1, 5, 10);
cvll    = nan(size(lam_L2));
fprintf('CV over lambda (held-out log-lik/obs, kernel peak):\n');
for il = 1:numel(lam_L2)
    cvll(il) = pen_grouped_cv(X, yev, lam_L2(il)*Pk, rowfold, cv_folds);
    bt = penalized_logit(X, yev, lam_L2(il)*Pk);
    fprintf('  lambda=%9.3g  cvLL=%.6f  kpeak=%.4f\n', lam_L2(il), cvll(il), max(abs(B*bt(ker))));
end
[~, ib] = max(cvll);  lambda = lam_L2(ib);        % best CV deviance
[bhat_p, Vp] = penalized_logit(X, yev, lambda * Pk);

% ── CV curve: held-out log-lik vs lambda (how flat/peaked the choice is) ──────
cvfh = figure('Color','w','Position',[80 80 640 470]); axcv = axes(cvfh); hold(axcv,'on');
plot(axcv, lam_L2, cvll, '-o', 'Color',[0.11 0.42 0.68], ...
    'MarkerFaceColor',[0.11 0.42 0.68], 'LineWidth',1.8, 'MarkerSize',6);
plot(axcv, lambda, cvll(ib), 'p', 'MarkerSize',18, ...
    'MarkerFaceColor',[0.85 0.37 0.12], 'MarkerEdgeColor','k');
xline(axcv, lambda, ':', 'Color',[0.85 0.37 0.12], 'LineWidth',1.2);
set(axcv,'XScale','log','Box','off','TickDir','out','Layer','top','FontSize',14);
xlabel(axcv,'smoothness  \lambda'); ylabel(axcv,'held-out log-likelihood / obs');
title(axcv, sprintf('grouped %d-fold CV   (selected \\lambda = %.3g)', cv_folds, lambda), ...
    'FontWeight','normal');
if save_out
    cvn = fullfile(fileparts(mfilename('fullpath')), sprintf('kernel_fit_cvcurve%s%s.png', ...
        tern(~include_baseline,'_nobase',''), tern(~include_nuisance,'_nonuis','')));
    exportgraphics(cvfh, cvn, 'Resolution', 200);
    [~, cvb] = fileparts(cvn);  exporter(cvfh, paths, [cvb '.pdf']);
    fprintf('Saved CV curve: %s\n', cvn);
end

% (b) unpenalised: ML fit in the bump basis (the projection is the only regulariser)
fprintf('Unpenalised ML fit in the %d-bump basis (lambda = 0).\n', K);
[bhat_u, Vu] = penalized_logit(X, yev, 0*Pk);

% reconstruct both kernels β(τ)=B·b with their CIs
bk_p     = bhat_p(ker);  kernel_p = B * bk_p;  kse_p = sqrt(sum((B*Vp(ker,ker)).*B, 2));
bk_u     = bhat_u(ker);  kernel_u = B * bk_u;  kse_u = sqrt(sum((B*Vu(ker,ker)).*B, 2));
fprintf('Penalised:   lambda=%.3g, ||b||=%.3f, peak=%.3f\n', lambda, norm(bk_p), max(abs(kernel_p)));
fprintf('Unpenalised: lambda=0,    ||b||=%.3f, peak=%.3f\n', norm(bk_u), max(abs(kernel_u)));

% colours for the two fits (shared by panels 2 & 3)
c_pen = [0.11 0.42 0.68];   % penalised  (blue)
c_un  = [0.85 0.37 0.10];   % unpenalised (orange)

% ── colours + bump centres (shared across the three panels) ──────────────────
cmapB  = cbrewer2('Spectral', K);
c_lag  = linspace(min(lag_fr), max(lag_fr), K);   % bump centres (frames)
c_trel = -c_lag / fps;                            % centre time (s)
[~, bord] = sort(c_trel);                         % past -> future

% ════════════════════════════════════════════════════════════════════════
% figure:  [ Xk | y ]   -fit->   [ b ]   -β=B·b->   [ β(τ) ]
% ════════════════════════════════════════════════════════════════════════
fh = figure('Color','w','Position',[50 60 1460 470]);
tl = tiledlayout(1, 3, 'TileSpacing','loose', 'Padding','compact'); %#ok<NASGU>
redmap = cbrewer2('Reds', []);
xl = [-kernel_past_s kernel_future_s];

% ── (1) input: a slice of Xk + the outcome y ─────────────────────────────────
ev = find(yev==1);
if ~isempty(ev), e = ev(ceil(numel(ev)/2)); rows_show = max(1,e-34):min(r,e+4);
else,           rows_show = 1:min(40,r); end
ns = numel(rows_show);
axA = nexttile(1); hold(axA,'on');
imagesc(axA, 1:K, 1:ns, Xk(rows_show, bord), 'AlphaData', ~isnan(Xk(rows_show,bord)));
set(axA,'YDir','reverse','Color',[0.86 0.86 0.86]);
colormap(axA, redmap);
qx = max(1, prctile(abs(Xk(:)), 99)); caxis(axA,[-qx qx]);
% outcome y as a column of squares just right of Xk
yshow = yev(rows_show);
ycol  = repmat([0.90 0.90 0.87], ns, 1); ycol(yshow==1,:) = repmat([0.94 0.63 0.15], sum(yshow==1), 1);
scatter(axA, repmat(K+1.2, ns, 1), 1:ns, 100, ycol, 'filled', 's', 'MarkerEdgeColor',[0.5 0.5 0.5]);
text(axA, K+1.2, -0.5, 'Break?', 'HorizontalAlignment','center', 'FontSize',16, 'VerticalAlignment', 'bottom');
xlim(axA,[0.5 K+2]); ylim(axA,[0.5 ns+0.5]);
set(axA,'XTick',[],'YTick',[]);
xlabel(axA,'basis columns'); ylabel(axA,'at-risk frames');
apply_generic(axA,'font_size',24);
lbl = cell(1,K);
for p = 1:K
    col = cmapB(bord(p),:);
    lbl{p} = sprintf('\\color[rgb]{%.3f,%.3f,%.3f}%.1f', col(1), col(2), col(3), c_trel(bord(p)));
end
set(axA,'XTick',1:K,'XTickLabel',lbl,'TickLabelInterpreter','tex','YTick',[]);

% ── (2) fitted coefficients b — penalised vs unpenalised (grouped bars) ──────
axB = nexttile(2); hold(axB,'on');
sp  = mean(diff(sort(c_trel)));
w   = 0.40*sp;  off = 0.22*sp;
for j = 1:K
    bar(axB, c_trel(j)-off, bk_p(j), w, 'FaceColor', c_pen, 'EdgeColor','none');
    bar(axB, c_trel(j)+off, bk_u(j), w, 'FaceColor', c_un,  'EdgeColor','none');
end
plot(axB, xl, [0 0], 'k-', 'LineWidth', 0.8);
hp = patch(axB, nan, nan, c_pen, 'EdgeColor','none');   % legend proxies
hu = patch(axB, nan, nan, c_un,  'EdgeColor','none');
xlim(axB, xl);
xlabel(axB,'bump time (s)'); ylabel(axB,'coefficient b_k');
apply_generic(axB,'font_size',24);

% ── (3) reconstructed kernel β(τ) = B·b — penalised vs unpenalised, with CIs ──
[t_ax, ord] = sort(t_rel);
axC = nexttile(3); hold(axC,'on');
kup_u = kernel_u + 1.96*kse_u;  klo_u = kernel_u - 1.96*kse_u;
kup_p = kernel_p + 1.96*kse_p;  klo_p = kernel_p - 1.96*kse_p;
patch(axC, [t_ax; flipud(t_ax)], [kup_u(ord); flipud(klo_u(ord))], c_un, ...
    'FaceAlpha', 0.12, 'EdgeColor','none');
patch(axC, [t_ax; flipud(t_ax)], [kup_p(ord); flipud(klo_p(ord))], c_pen, ...
    'FaceAlpha', 0.15, 'EdgeColor','none');
plot(axC, xl, [0 0], 'k-', 'LineWidth', 0.8);
hcu = plot(axC, t_ax, kernel_u(ord), '-', 'Color', c_un,  'LineWidth', 2.0);
hcp = plot(axC, t_ax, kernel_p(ord), '-', 'Color', c_pen, 'LineWidth', 2.6);
legend(axC, [hcp hcu], {'penalised','unpenalised'}, 'Box','off','Location','northoutside', 'FontSize',  24);
xlim(axC, xl);
ylim([-0.02 0.02])
xlabel(axC,'time relative to un-freeze (s)'); ylabel(axC,'\beta(\tau)');
apply_generic(axC,'font_size',24);

% ── flow arrows between the panels (in the gaps, labels above the line) ──────
drawnow;
gapA = mean([axA.Position(1)+axA.Position(3), axB.Position(1)]);   % gap 1 centre
gapB = mean([axB.Position(1)+axB.Position(3), axC.Position(1)]);   % gap 2 centre
ay   = axA.Position(2) + 0.55*axA.Position(4);


% ── save ─────────────────────────────────────────────────────────────────────
set([axA axB axC],'Layer','top');
axA.Toolbar.Visible='off'; axB.Toolbar.Visible='off'; axC.Toolbar.Visible='off';
if save_out
    outdir = fileparts(mfilename('fullpath'));
    fn = fullfile(outdir, sprintf('kernel_fit_compare_%s%s%s%s.png', grid_anchor, ...
        tern(mask_preonset,'_masked',''), tern(~include_baseline,'_nobase',''), tern(~include_nuisance,'_nonuis','')));
    exportgraphics(fh, fn, 'Resolution', 240);
    [~, base] = fileparts(fn);
    exporter(fh, paths, [base '.pdf']);   % vector PDF -> paths.fig (repo convention)
    fprintf('\nSaved: %s\n       %s\n', fn, fullfile(paths.fig, [base '.pdf']));
end

% ════════════════════════════════════════════════════════════════════════
% local functions (copied verbatim from hazard_kernel_sparse.m)
% ════════════════════════════════════════════════════════════════════════
function [lambda, cvll] = select_lambda_l2(X, y, Pk, lam_grid, rowfold, K)
    cvll_lam = nan(size(lam_grid));  cvse_lam = nan(size(lam_grid));
    for il = 1:numel(lam_grid)
        [cvll_lam(il), llf] = pen_grouped_cv(X, y, lam_grid(il)*Pk, rowfold, K);
        cvse_lam(il) = std(llf, 'omitnan') / sqrt(sum(~isnan(llf)));
    end
    [~, ibest] = max(cvll_lam);
    thr   = cvll_lam(ibest) - cvse_lam(ibest);
    i1se  = find(cvll_lam >= thr, 1, 'last');
    lambda = lam_grid(i1se);
    cvll   = cvll_lam(i1se);
end

function D = diffmat(n, ord)
    D = eye(n);
    for k = 1:ord, D = diff(D); end
end

function [b, V] = penalized_logit(X, y, P)
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

function [ll, ll_folds] = pen_grouped_cv(X, y, P, rowfold, K)
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

function s = tern(c,a,b), if c, s=a; else, s=b; end, end
