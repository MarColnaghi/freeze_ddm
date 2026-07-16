clearvars; clc;
% ════════════════════════════════════════════════════════════════════════
% kernel_fit_nbasis — VISUAL AID  (basis-invariance check)
%
% Refit the SAME hazard model with a different number of raised-cosine bumps
% (K = 6, 12, 24) and overlay the recovered kernel β(τ). The point:
%
%   • the fitted coefficients b live in the basis and change a lot with K,
%   • but the KERNEL β(τ) = B·b lives on the lag axis and is basis-invariant
%     — provided K is large enough to represent the shape.
%
% One panel overlays K = 6, 12, 24 for a single fit mode set by `smoothness`:
%   smoothness = true  → CV-tuned L2 penalty (curves sit almost on top of each other),
%   smoothness = false → unsmoothed ML in the bump basis (more wiggle as K grows).
% ════════════════════════════════════════════════════════════════════════

% ── parameters (mirror kernel_fit.m) ─────────────────────────────────────────
fps             = 60;
kernel_past_s   = 4.0;
kernel_future_s = 2;
dt_frames       = 6;
trunc_point     = 30;   entry_fr = trunc_point;
contact_threshold = 70;
grid_anchor     = 'offset';
mask_preonset   = true;

nb_list         = [6 12 24];   % kernel bases to compare
bump_width      = 1.25;
nb_baseline     = 6;           % baseline-hazard bumps over time-in-freeze (nuisance)
cv_folds        = 5;           % grouped CV folds to pick the smoothness penalty
smoothness      = false;       % false → UNSMOOTHED (unpenalised ML fit in the bump basis)
include_baseline = false;      % false → drop the time-in-freeze baseline hazard block
include_nuisance = false;      % false → drop the moving_flies + loom-speed dummies
save_out        = true;
rng(11123);

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

% ── build the FULL person-period design S, outcome yev, nuisances (ONCE) ─────
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

% ── nuisances + shared CV folds (same for every K, so fits are comparable) ───
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
fprintf('Design cols: intercept + K kernel + %d nuisance (baseline=%d, mov/slo=%d)\n', ...
    size(Nuis,2), include_baseline, include_nuisance);

fold    = mod(randperm(max(gid)), cv_folds) + 1;
rowfold = fold(gid);
lam_L2  = logspace(-1, 5, 10);

% ── refit for each K in the selected mode, store one kernel ──────────────────
nK      = numel(nb_list);
Kk      = nan(nLag, nK);   Kse = nan(nLag, nK);   % kernel + SE
lam_sel = nan(1, nK);
for ik = 1:nK
    nb = nb_list(ik);
    B  = raised_cosine_basis(lag_fr, nb, bump_width);   % nLag × K
    K  = size(B,2);
    Xk = S * B;
    X  = [ones(r,1), Xk, Nuis];
    ker = 1 + (1:K);
    Dk  = diffmat(K, 2);
    Pk  = zeros(size(X,2)); Pk(ker,ker) = Dk' * Dk;

    if smoothness
        cvll = nan(size(lam_L2));                        % pick lambda by grouped CV
        for il = 1:numel(lam_L2)
            cvll(il) = pen_grouped_cv(X, yev, lam_L2(il)*Pk, rowfold, cv_folds);
        end
        [~, ib] = max(cvll);  lam_sel(ik) = lam_L2(ib);
    else
        lam_sel(ik) = 0;                                 % UNSMOOTHED: unpenalised ML
    end
    [bf, Vf]  = penalized_logit(X, yev, lam_sel(ik) * Pk);
    Kk(:,ik)  = B * bf(ker);
    Kse(:,ik) = sqrt(sum((B*Vf(ker,ker)).*B, 2));
    fprintf('K=%2d: lambda=%9.3g  peak=%.4f\n', nb, lam_sel(ik), max(abs(Kk(:,ik))));
end

% ── how basis-invariant is the shape? corr each kernel with the K=24 fit ─────
ref = nK;   % reference = richest basis
fprintf('\nShape agreement vs K=%d (Pearson r over all lags):\n', nb_list(ref));
for ik = 1:nK
    fprintf('  K=%2d:  r=%.4f\n', nb_list(ik), corr(Kk(:,ik), Kk(:,ref)));
end

% ════════════════════════════════════════════════════════════════════════
% figure: β(τ) overlaid across K — penalised (L) vs unpenalised (R)
% ════════════════════════════════════════════════════════════════════════
cmapK = [0.20 0.45 0.75;    % K=6  blue
         0.25 0.62 0.32;    % K=12 green
         0.85 0.37 0.12];   % K=24 orange
[t_ax, ord] = sort(t_rel);
xl = [-kernel_past_s kernel_future_s];

fh = figure('Color','w','Position',[50 60 500 500]);
tl = tiledlayout(1, 1, 'TileSpacing','compact', 'Padding','compact'); %#ok<NASGU>


% ── β(τ) overlaid across K (selected fit mode) ───────────────────────────────
axU = nexttile(1); hold(axU,'on');
hU = gobjects(1,nK);
plot(axU, xl, [0 0], 'k-', 'LineWidth', 0.8);
for ik = 1:nK
    hU(ik) = plot(axU, t_ax, Kk(ord,ik), '-', 'Color', cmapK(ik,:), 'LineWidth', 2.0);
end
xlim(axU, xl);
xlabel(axU,'time relative to un-freeze (s)'); ylabel(axU,'\beta(\tau)');
legend(axU, hU, compose('K = %d', nb_list), 'Box','off','Location','best', 'FontSize', 20);
apply_generic(axU,'font_size', 24); axU.Toolbar.Visible='off';

% ── save ─────────────────────────────────────────────────────────────────────
if save_out
    outdir = fileparts(mfilename('fullpath'));
    fn = fullfile(outdir, sprintf('kernel_fit_nbasis_%s%s%s%s%s.png', grid_anchor, ...
        tern(mask_preonset,'_masked',''), tern(~smoothness,'_nosmooth',''), ...
        tern(~include_baseline,'_nobase',''), tern(~include_nuisance,'_nonuis','')));
    exportgraphics(fh, fn, 'Resolution', 220);   % PNG preview
    [~, base] = fileparts(fn);
    exporter(fh, paths, [base '.pdf']);   % vector PDF -> paths.fig (repo convention)
    fprintf('\nSaved: %s\n       %s\n', fn, fullfile(paths.fig, [base '.pdf']));
end

% ════════════════════════════════════════════════════════════════════════
% local functions (copied verbatim from kernel_fit.m)
% ════════════════════════════════════════════════════════════════════════
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
