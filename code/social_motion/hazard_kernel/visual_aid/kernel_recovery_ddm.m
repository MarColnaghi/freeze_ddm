clearvars; clc;
% ════════════════════════════════════════════════════════════════════════
% kernel_recovery_ddm — VISUAL AID / VALIDATION
%
% Ground-truth recovery. Instead of the REAL un-freeze times, we SIMULATE them
% from a bounded (leaky) accumulator with FIXED parameters, driven by the REAL
% social motion each fly saw:
%
%     x(k+1) = x(k)(1 - λΔt) + (μ0 + β·sm(k))Δt + σ√Δt·randn,   λ = 1/τ_gen
%     un-freeze = first k with x ≥ θ         (sim_leaky_accumulator)
%
% The generative leak gives a KNOWN hazard kernel: β_gt(τ) ∝ exp(-τ/τ_gen) on
% causal lags (τ_gen = Inf → flat/perfect DDM), and ≈0 on acausal lags. We then
% recover β(τ) with the SAME machinery as kernel_fit_single.m (raised-cosine
% basis + penalised/unpenalised logistic) and overlay the ground truth. If the
% recovered kernel matches the exponential shape (and stays ~0 acausally), the
% estimator is trustworthy.
%
% Uses the existing generator (generate_accumulator_dataset + sim_leaky_accumulator).
% ════════════════════════════════════════════════════════════════════════

% ── paths to the simulator + generator ───────────────────────────────────────
here     = fileparts(mfilename('fullpath'));                 % .../visual_aid
hk       = fileparts(here);                                  % .../hazard_kernel
code_dir = fileparts(fileparts(fileparts(hk)));              % .../code
addpath(fullfile(hk, 'recovery'));                           % generate_accumulator_dataset
addpath(fullfile(code_dir, 'sims', 'simulators'));           % sim_leaky_accumulator

% ── recovery params (mirror kernel_fit_single.m) ─────────────────────────────
fps             = 60;
kernel_past_s   = 4.0;
kernel_future_s = 2;      % acausal control — should recover ≈ 0
dt_frames       = 6;
trunc_point     = 30;   entry_fr = trunc_point;
grid_anchor     = 'offset';
mask_preonset   = true;

nb_basis        = 12;
bump_width      = 1.25;
nb_baseline     = 60;
cv_folds        = 20;
smoothness      = true;    % recover a smooth kernel (true) or unsmoothed ML (false)
include_baseline = false;
include_nuisance = false;
save_out        = true;
rng(722);

% ── GROUND-TRUTH generative DDM (fixed parameters) ───────────────────────────
tau_gen   = Inf;           % generative leak time-constant (s); Inf → perfect DDM (flat kernel)
P.dt      = 1/fps;
P.mu0     = 1.6;             % baseline drive
P.beta    = 5;             % social gain
P.sigma   = 1.0;           % diffusion SD
P.theta   = 6;             % bound
P.n_aug   = 5;             % synthetic trajectories per real bout (more events)
P.maxT_fr = round(10.5*fps);
P.trunc_point = trunc_point;
P.seed    = 1;
P.lambda  = 0; if ~isinf(tau_gen), P.lambda = 1/tau_gen; end

% ── load real bouts + motion cache (template for the generator) ──────────────
threshold_imm = 2; threshold_mob = 2; threshold_pc = 4;
id_code = sprintf('imm%d_mob%d_pc%d', threshold_imm, threshold_mob, threshold_pc);
paths        = path_generator('folder','social_motion/hazard_kernel/visual_aid','bouts_id',id_code,'imfirst',false);
bouts        = importdata(fullfile(paths.dataset,  'bouts.mat'));
motion_cache = importdata(fullfile(paths.cache_path,'motion_cache.mat'));
thresholds = define_thresholds('le_window', struct('le_window_sl',[5 55],'le_window_fl',[5 55]));
bouts = bouts_formatting(bouts, thresholds);
bl_t  = data_parser_new(bouts, 'type','immobility','period','loom','window','le','nloom',2:20,'min_dur',trunc_point);
nB = height(bl_t);
fprintf('Template: %d real freezes.\n', nB);

% ── per-bout sm drive for the generator (onset forward, capped, /10, NaN->0) ─
sm_in = cell(nB,1);
for b = 1:nB
    sm  = motion_cache(bl_t.fly(b)); sm = sm(:) ./ 10;
    on  = bl_t.onsets(b);
    cap = min(P.maxT_fr, numel(sm) - on + 1);
    sm_in{b} = sm(on : on + cap - 1);
    sm_in{b}(isnan(sm_in{b})) = 0;
end

% ── SIMULATE freezes from the DDM (fixed params, real sm drive) ──────────────
[bl, gen_info] = generate_accumulator_dataset(bl_t, sm_in, P);
fprintf('Simulated: %d synthetic freezes (%.1f%% censored), median dur %.2f s\n', ...
    height(bl), 100*gen_info.cens_frac, gen_info.median_fr/fps);

% ── global z-scoring of sm (for the estimator's design) ──────────────────────
flies = unique(bl_t.fly); nn=0; s1=0; s2=0;
for i = 1:numel(flies)
    v = motion_cache(flies(i)); v = v(:); v = v(~isnan(v));
    nn=nn+numel(v); s1=s1+sum(v); s2=s2+sum(v.^2);
end
sm_mu = s1/nn; sm_sd = sqrt(s2/nn - sm_mu^2);

% ── lag grid ─────────────────────────────────────────────────────────────────
back_fr = round(kernel_past_s*fps);  fwd_fr = round(kernel_future_s*fps);
lag_fr  = (-fwd_fr:back_fr)';  nLag = numel(lag_fr);  lag_s = lag_fr/fps;
t_rel   = -lag_s;  caus_mask = lag_fr >= 0;

% ── build the person-period design S from the SYNTHETIC bl + real sm ─────────
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
B  = raised_cosine_basis(lag_fr, nb_basis, bump_width);
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

% ── fit exactly as kernel_fit_single.m (penalised or unsmoothed) ─────────────
X   = [ones(r,1), Xk, Nuis];
ker = 1 + (1:K);
Dk  = diffmat(K, 2);
Pk  = zeros(size(X,2)); Pk(ker,ker) = Dk' * Dk;
if smoothness
    fold    = mod(randperm(max(gid)), cv_folds) + 1;
    rowfold = fold(gid);
    lam_L2  = logspace(-1, 5, 10);
    cvll    = nan(size(lam_L2));
    for il = 1:numel(lam_L2)
        cvll(il) = pen_grouped_cv(X, yev, lam_L2(il)*Pk, rowfold, cv_folds);
    end
    [~, ib] = max(cvll);  lambda = lam_L2(ib);
else
    lambda = 0;
end
[bhat, V] = penalized_logit(X, yev, lambda * Pk);
bk     = bhat(ker);
kernel = B * bk;
kse    = sqrt(sum((B * V(ker,ker)) .* B, 2));

% ── ground-truth kernel (generative): exp(-τ/τ_gen) causal, 0 acausal ────────
gt = zeros(nLag,1);
if isinf(tau_gen), gt(caus_mask) = 1; else, gt(caus_mask) = exp(-lag_s(caus_mask)/tau_gen); end
c   = caus_mask;                                  % least-squares scale to the recovered kernel
a_s = (gt(c)' * kernel(c)) / max(gt(c)' * gt(c), eps);
gt_scaled = a_s * gt;

% ── recovered timescale (exp fit to the causal kernel) vs generative ─────────
tau_hat = fit_leak_tau(kernel, lag_s, caus_mask, kse);
fprintf('lambda=%.3g  |  generative tau=%.2fs  recovered tau_hat=%.2fs  |  peak=%.3f  acausal max|beta|=%.3f\n', ...
    lambda, tau_gen, tau_hat, max(abs(kernel(caus_mask))), max(abs(kernel(~caus_mask))));

% ── colours + bump centres (as kernel_fit_single.m) ──────────────────────────
cmapB  = cbrewer2('Spectral', K);
c_lag  = linspace(min(lag_fr), max(lag_fr), K);
c_trel = -c_lag / fps;
[~, bord] = sort(c_trel);

% ════════════════════════════════════════════════════════════════════════
% figure:  [ Xk | y ]   -fit->   [ b ]   -β=B·b->   [ β(τ): recovered vs truth ]
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
yshow = yev(rows_show);
ycol  = repmat([0.90 0.90 0.87], ns, 1); ycol(yshow==1,:) = repmat([0.94 0.63 0.15], sum(yshow==1), 1);
scatter(axA, repmat(K+1.2, ns, 1), 1:ns, 100, ycol, 'filled', 's', 'MarkerEdgeColor',[0.5 0.5 0.5]);
text(axA, K+1.2, -0.5, 'Break?', 'HorizontalAlignment','center', 'FontSize',16, 'VerticalAlignment','bottom');
xlim(axA,[0.5 K+2]); ylim(axA,[0.5 ns+0.5]);
set(axA,'XTick',[],'YTick',[]);
xlabel(axA,'basis columns'); ylabel(axA,'at-risk frames');
apply_generic(axA,'font_size',20);
lbl = cell(1,K);
for p = 1:K
    col = cmapB(bord(p),:);
    lbl{p} = sprintf('\\color[rgb]{%.3f,%.3f,%.3f}%.1f', col(1), col(2), col(3), c_trel(bord(p)));
end
set(axA,'XTick',1:K,'XTickLabel',lbl,'TickLabelInterpreter','tex','YTick',[]);

% ── (2) fitted coefficients b (one bar per bump, coloured by basis) ──────────
axB = nexttile(2); hold(axB,'on');
sp = mean(diff(sort(c_trel)));
for j = 1:K
    bar(axB, c_trel(j), bk(j), 0.85*sp, 'FaceColor', cmapB(j,:), 'EdgeColor','none');
end
plot(axB, xl, [0 0], 'k-', 'LineWidth', 0.8);
xlim(axB, [xl(1) - 1; xl(2) + 1 ]);
xlabel(axB,'bump time (s)'); ylabel(axB,'coefficient b_k');
apply_generic(axB,'font_size',20);

% ── (3) recovered β(τ) + CI vs the generative ground-truth kernel ────────────
[t_ax, ord] = sort(t_rel);
axC = nexttile(3); hold(axC,'on');
ku = kernel + 1.96*kse; klw = kernel - 1.96*kse;
patch(axC, [t_ax; flipud(t_ax)], [ku(ord); flipud(klw(ord))], [0.2 0.2 0.2], ...
    'FaceAlpha', 0.15, 'EdgeColor','none');
plot(axC, xl, [0 0], 'k-', 'LineWidth', 0.8);
hrec = plot(axC, t_ax, kernel(ord),    '-',  'Color',[0.10 0.10 0.10], 'LineWidth', 2.6);
xlim(axC, xl);
xlabel(axC,'time relative to un-freeze (s)'); ylabel(axC,'\beta(\tau)');
apply_generic(axC,'font_size',20);

set([axA axB axC],'Layer','top');
axA.Toolbar.Visible='off'; axB.Toolbar.Visible='off'; axC.Toolbar.Visible='off';

% ── save ─────────────────────────────────────────────────────────────────────
if save_out
    outdir = fileparts(mfilename('fullpath'));
    gtag = ternstr(isinf(tau_gen), 'flat', sprintf('tau%.2f', tau_gen));
    fn = fullfile(outdir, sprintf('kernel_recovery_ddm_%s_%s%s.png', grid_anchor, gtag, tern(~smoothness,'_nosmooth','')));
    exportgraphics(fh, fn, 'Resolution', 240);
    [~, base] = fileparts(fn);
    exporter(fh, paths, [base '.pdf']);   % vector PDF -> paths.fig (repo convention)
    fprintf('\nSaved: %s\n       %s\n', fn, fullfile(paths.fig, [base '.pdf']));
end

% ════════════════════════════════════════════════════════════════════════
% local functions (copied from kernel_fit_single.m / recovery)
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

function tau_hat = fit_leak_tau(kernel, lag_s, caus, kernel_se)
% Weighted exponential fit β(τ) ≈ A·exp(-τ/τ)+c to the causal kernel.
    tau_s = lag_s(caus);  kc = kernel(caus);
    if nargin >= 4 && ~isempty(kernel_se) && any(kernel_se(caus) > 0)
        w = 1 ./ max(kernel_se(caus), eps).^2;
    else
        w = ones(size(kc));
    end
    expmod = @(p,t) p(1).*exp(-t./exp(p(2))) + p(3);
    sse    = @(p) sum(w .* (kc - expmod(p,tau_s)).^2);
    p0     = [max(kc)-min(kc), log(0.5), min(kc)];
    opt    = optimset('Display','off','MaxFunEvals',1e4,'MaxIter',1e4);
    try, pf = fminsearch(sse, p0, opt); tau_hat = exp(pf(2)); catch, tau_hat = NaN; end
end

function s = tern(c,a,b), if c, s=a; else, s=b; end, end
function s = ternstr(c,a,b), if c, s=a; else, s=b; end, end
