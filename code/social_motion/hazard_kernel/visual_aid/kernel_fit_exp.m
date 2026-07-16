clearvars; clc;
% ════════════════════════════════════════════════════════════════════════
% kernel_fit_exp — VISUAL AID  (exponential-bank companion to kernel_fit.m)
%
% Same loop — [ Xₖ | y ] --fit--> b --β=B·b--> β(τ) — but with the EXPONENTIAL
% BANK and an L1 (sparse) fit, i.e. model (B) of hazard_kernel_sparse.m:
%
%   • B = exp_basis: decaying exponentials exp(-τ/τ_k), log-spaced time-constants.
%   • L1 penalty ⇒ most coefficients are driven to EXACTLY zero; the survivors ARE
%     the selected timescales (the leak). MIDDLE panel: the sparse coefficients b
%     — a few nonzero bars at the retained τ_k.
%   • RIGHT panel: β(τ) = B·b is the sum of just those few weighted decays (+ CI).
%
% Fit uses ALL at-risk intervals; the LEFT panel shows a slice. Toggle
% include_baseline to add/remove the time-in-freeze baseline hazard.
% ════════════════════════════════════════════════════════════════════════

% ── parameters (mirror hazard_kernel_sparse.m model B) ───────────────────────
fps             = 60;
kernel_past_s   = 4.0;
kernel_future_s = 0;      % model B uses no acausal block (pure causal exponentials)
dt_frames       = 6;
trunc_point     = 30;   entry_fr = trunc_point;
contact_threshold = 70;
grid_anchor     = 'offset';
mask_preonset   = false;  % match the sibling kernel_fit.m; flip to true for within-freeze-only

n_exp           = 10;     % exponential atoms
tau_min_s       = 0.05;   % fastest time-constant (s)
tau_max_s       = 8.0;    % slowest time-constant (s)
nb_acausal      = 0;      % acausal control bumps (0 with kernel_future_s = 0)
bump_width      = 1.25;
nb_baseline     = 6;      % baseline-hazard bumps over time-in-freeze
include_baseline = false;  % false → drop the time-in-freeze baseline hazard block
include_nuisance = false;  % false → drop the moving_flies + loom-speed dummies
cv_folds        = 20;
sel_rule        = 'argmax';  % lambda pick: 'argmax' (best CV) | '1se' (sparser: largest lambda within 1 SE)
save_out        = true;
rng(12320);

% ── load + filter (identical to the real pipeline) ───────────────────────────
threshold_imm = 2; threshold_mob = 2; threshold_pc = 4;
id_code = sprintf('imm%d_mob%d_pc%d', threshold_imm, threshold_mob, threshold_pc);
paths        = path_generator('folder','social_motion/hazard_kernel/visual_aid','bouts_id',id_code,'imfirst',false);
bouts        = importdata(fullfile(paths.dataset,  'bouts.mat'));
motion_cache = importdata(fullfile(paths.cache_path,'motion_cache.mat'));
thresholds = define_thresholds('le_window', struct('le_window_sl',[5 55],'le_window_fl',[5 55]));
bouts = bouts_formatting(bouts, thresholds);

bl = data_parser_new(bouts, 'type','immobility','period','loom','window','le','nloom',2:20,'min_dur',trunc_point);
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

% ── exponential bank, projection, nuisances ──────────────────────────────────
tau_grid_s = logspace(log10(tau_min_s), log10(tau_max_s), n_exp);
[B, tau_grid, is_exp] = exp_basis(lag_fr, tau_grid_s, fps, nb_acausal, bump_width);
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
fprintf('Design cols: intercept + %d atoms + %d nuisance (baseline=%d, mov/slo=%d)\n', ...
    K, size(Nuis,2), include_baseline, include_nuisance);

% ── L1 fit with a grouped-CV λ-path (argmax CV), fast covariance-update solver ─
X       = [ones(r,1), Xk, Nuis];
ker     = 1 + (1:K);
penmask = false(1, size(X,2)); penmask(ker) = true;

b_nuis  = penalized_logit(X(:,~penmask), yev, 0);
mu0     = 1 ./ (1 + exp(-min(max(X(:,~penmask)*b_nuis,-30),30)));
lam_max = max(abs(X(:,penmask)' * (yev - mu0)));
lam_grid = lam_max * logspace(0, -3, 12);

fold    = mod(randperm(max(gid)), cv_folds) + 1;
rowfold = fold(gid);
llf = nan(cv_folds, numel(lam_grid));
for k = 1:cv_folds
    fprintf('  CV fold %d/%d...\n', k, cv_folds);
    te = rowfold==k; tr = ~te; bw = [];
    for il = 1:numel(lam_grid)
        bw = penalized_logit_l1(X(tr,:), yev(tr), lam_grid(il), penmask, bw);
        eta = min(max(X(te,:)*bw, -30), 30);
        mu  = min(max(1./(1+exp(-eta)), 1e-9), 1-1e-9);
        llf(k,il) = mean(yev(te).*log(mu) + (1-yev(te)).*log(1-mu));
    end
end
cvll_lam = mean(llf, 1, 'omitnan');
[~, isel] = max(cvll_lam);  lambda = lam_grid(isel);
% (a) penalised (sparse L1, CV lambda) + debiased CI on the selected atoms
bhat_p   = penalized_logit_l1(X, yev, lambda, penmask, []);
bk_p     = bhat_p(ker);
kernel_p = B * bk_p;
active   = find(abs(bk_p) > 1e-8);
if isempty(active)
    kse_p = zeros(size(kernel_p));  sel_tau = [];
else
    Xdb = [ones(r,1), S*B(:,active), Nuis];
    [~, Vdb] = penalized_logit(Xdb, yev, 1e-6*eye(size(Xdb,2)));
    kcol  = 1 + (1:numel(active));
    kse_p = sqrt(sum((B(:,active) * Vdb(kcol,kcol)) .* B(:,active), 2));
    sel_tau = tau_grid(active(is_exp(active)));
end

% (b) unpenalised ML in the exp bank (lambda = 0; the projection is the only regulariser)
[bhat_u, Vu] = penalized_logit(X, yev, 1e-6*eye(size(X,2)));
bk_u     = bhat_u(ker);
kernel_u = B * bk_u;
kse_u    = sqrt(sum((B*Vu(ker,ker)).*B, 2));

fprintf('Penalised L1: lambda=%.3g, %d/%d atoms selected (peak %.3f)\n', lambda, numel(active), K, max(abs(kernel_p)));
fprintf('Unpenalised : peak %.3f\n', max(abs(kernel_u)));
if ~isempty(sel_tau)
    fprintf('Selected tau (s): %s\n', strjoin(arrayfun(@(t)sprintf('%.2f',t), sel_tau,'uni',0), ', '));
end

% ── colours: exp-atom Spectral (for τ ticks) + the two fits (bars/curves) ────
cmapE = cbrewer2('Spectral', n_exp);
c_pen = [0.11 0.42 0.68];   % penalised L1 (blue)
c_un  = [0.85 0.37 0.10];   % unpenalised (orange)

% ════════════════════════════════════════════════════════════════════════
% figure:  [ Xk | y ]   -fit->   [ b (sparse) ]   -β=B·b->   [ β(τ) ]
% ════════════════════════════════════════════════════════════════════════
fh = figure('Color','w','Position',[50 60 1460 470]);
tl = tiledlayout(1, 3, 'TileSpacing','loose', 'Padding','compact'); %#ok<NASGU>
redmap = cbrewer2('Reds', []);
xl = [-kernel_past_s max(0, kernel_future_s)];

% ── (1) input: a slice of Xk + the outcome y ─────────────────────────────────
ev = find(yev==1);
if ~isempty(ev), e = ev(ceil(numel(ev)/2)); rows_show = max(1,e-34):min(r,e+4);
else,           rows_show = 1:min(40,r); end
ns = numel(rows_show);
axA = nexttile(1); hold(axA,'on');
imagesc(axA, 1:K, 1:ns, Xk(rows_show,:), 'AlphaData', ~isnan(Xk(rows_show,:)));
set(axA,'YDir','reverse','Color',[0.86 0.86 0.86]);
colormap(axA, redmap);
qx = max(1, prctile(abs(Xk(:)), 99)); caxis(axA,[-qx qx]);
yshow = yev(rows_show);
ycol  = repmat([0.90 0.90 0.87], ns, 1); ycol(yshow==1,:) = repmat([0.94 0.63 0.15], sum(yshow==1), 1);
scatter(axA, repmat(K+1.2, ns, 1), 1:ns, 100, ycol, 'filled', 's', 'MarkerEdgeColor',[0.5 0.5 0.5]);
xlim(axA,[0.5 K+2]); ylim(axA,[0.5 ns+0.5]);
set(axA,'XTick',[],'YTick',[]);
xlabel(axA,'exp atoms  X_k'); ylabel(axA,'at-risk frames');
apply_generic(axA,'font_size',24);
lblA = repmat({''},1,K);
for j = 1:n_exp
    lblA{j} = sprintf('\\color[rgb]{%.3f,%.3f,%.3f}%.2f', cmapE(j,1),cmapE(j,2),cmapE(j,3), tau_grid(j));
end
set(axA,'XTick',1:K,'XTickLabel',lblA,'TickLabelInterpreter','tex','YTick',[]);
axA.XAxis.FontSize = 24;

% ── (2) coefficients b — penalised L1 (left axis) vs unpenalised (right axis) ─
% Two y-axes: the sparse penalised atoms and the large collinear unpenalised
% ones each get their own scale; symmetric limits keep both zeros aligned.
axB = nexttile(2); hold(axB,'on');
wb = 0.38; off = 0.2;
mp = max(abs(bk_p))*1.2 + eps;   mu = max(abs(bk_u))*1.2 + eps;
yyaxis(axB,'left');
for j = 1:K, bar(axB, j-off, bk_p(j), wb, 'FaceColor', c_pen, 'EdgeColor','none'); end
ylim(axB,[-mp mp]);  ylabel(axB,'penalised  b_k');
yyaxis(axB,'right');
for j = 1:K, bar(axB, j+off, bk_u(j), wb, 'FaceColor', c_un, 'EdgeColor','none'); end
ylim(axB,[-mu mu]);  ylabel(axB,'unpenalised  b_k');
plot(axB, [0.5 K+0.5], [0 0], '-', 'Color',[0.6 0.6 0.6], 'LineWidth', 0.6);
xlim(axB,[0.5 K+0.5]);
lbl = repmat({''},1,K);
for j = 1:n_exp, lbl{j} = sprintf('\\color[rgb]{%.3f,%.3f,%.3f}%.2f', cmapE(j,1),cmapE(j,2),cmapE(j,3), tau_grid(j)); end
xlabel(axB,'time-constant  \tau (s)');
yyaxis(axB,'left');  axB.YColor = c_pen;
yyaxis(axB,'right'); axB.YColor = c_un;
axB.YAxis(1).FontSize = 24;
axB.YAxis(2).FontSize = 24;
axB.YAxis(1).LineWidth = 2;
axB.YAxis(2).LineWidth = 2;
axB.XAxis.LineWidth = 2;

set(axB,'XTick',1:K,'XTickLabel',lbl,'TickLabelInterpreter','tex');
axB.XAxis.FontSize = 24;

% ── (3) reconstructed kernel β(τ) = B·b — penalised L1 vs unpenalised, with CIs ─
[t_ax, ord] = sort(t_rel);
axC = nexttile(3); hold(axC,'on');
kup_u = kernel_u + 1.96*kse_u; klo_u = kernel_u - 1.96*kse_u;
kup_p = kernel_p + 1.96*kse_p; klo_p = kernel_p - 1.96*kse_p;
patch(axC, [t_ax; flipud(t_ax)], [kup_u(ord); flipud(klo_u(ord))], c_un, ...
    'FaceAlpha', 0.12, 'EdgeColor','none');
patch(axC, [t_ax; flipud(t_ax)], [kup_p(ord); flipud(klo_p(ord))], c_pen, ...
    'FaceAlpha', 0.15, 'EdgeColor','none');
plot(axC, xl, [0 0], 'k-', 'LineWidth', 0.8);
hcu = plot(axC, t_ax, kernel_u(ord), '-', 'Color', c_un,  'LineWidth', 2.0);
hcp = plot(axC, t_ax, kernel_p(ord), '-', 'Color', c_pen, 'LineWidth', 2.6);
legend(axC, [hcp hcu], {'penalised L1','unpenalised'}, 'Box','off','Location','best', 'FontSize', 20 );
xlim(axC, xl);
xlabel(axC,'time relative to un-freeze (s)'); ylabel(axC,'\beta(\tau)');
apply_generic(axC,'font_size',24);

% ── save ─────────────────────────────────────────────────────────────────────
set([axA axB axC],'Layer','top');
axA.Toolbar.Visible='off'; axB.Toolbar.Visible='off'; axC.Toolbar.Visible='off';
if save_out
    outdir = fileparts(mfilename('fullpath'));
    bl_tag = ''; if ~include_baseline, bl_tag = '_nobase'; end
    nu_tag = ''; if ~include_nuisance, nu_tag = '_nonuis'; end
    fn = fullfile(outdir, sprintf('kernel_fit_exp_%s_%s%s%s.png', grid_anchor, sel_rule, bl_tag, nu_tag));
    exportgraphics(fh, fn, 'Resolution', 200);
    [~, base] = fileparts(fn);
    exporter(fh, paths, [base '.pdf']);   % vector PDF -> paths.fig (repo convention)
    fprintf('\nSaved: %s\n       %s\n', fn, fullfile(paths.fig, [base '.pdf']));
end

% ════════════════════════════════════════════════════════════════════════
% local functions (copied from hazard_kernel_sparse.m)
% ════════════════════════════════════════════════════════════════════════
function [B, tau_grid, is_exp] = exp_basis(lag_fr, tau_grid_s, fps, nb_acausal, bump_width)
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

function [b, support] = penalized_logit_l1(X, y, lambda, penmask, b0)
% Covariance-update (glmnet-style) L1 logistic — fast on collinear designs.
    [~, p] = size(X);
    if nargin < 5 || isempty(b0), b = zeros(p,1); else, b = b0(:); end
    penmask = logical(penmask(:));
    for outer = 1:100 %#ok<NASGU>
        eta = min(max(X*b, -30), 30);
        mu  = 1 ./ (1 + exp(-eta));
        w   = max(mu .* (1 - mu), 1e-6);
        z   = eta + (y - mu) ./ w;
        Xw  = X .* w;
        A   = Xw' * X;
        c   = Xw' * z;
        d   = max(diag(A), eps);
        b_old = b;  rc = c - A*b;
        for inner = 1:1000 %#ok<NASGU>
            bmax = 0;
            for j = 1:p
                rho = rc(j) + d(j)*b(j);
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
