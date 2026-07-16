function res = run_hazard_kernel(bl, motion_cache, opts)
% RUN_HAZARD_KERNEL  Discrete-time hazard kernel of freeze termination vs social motion.
%
%   res = run_hazard_kernel(bl, motion_cache, opts)
%
%   Analysis core extracted from hazard_kernel_sm.m so the real data and the
%   ground-truth simulations run through IDENTICAL methodology.
%
%   bl (table) must contain, one row per freeze bout:
%       fly          : fly id (key into motion_cache)
%       onsets       : absolute onset frame in the recording
%       dur_frames   : observed (possibly censored) freeze length, in frames
%       is_censored  : logical, true = right-censored (no observed un-freeze)
%       moving_flies : condition (used as a categorical control if >1 level)
%       sloom        : loom speed   (categorical control if >1 level)
%       bout_id      : unique id per bout (for grouped CV)
%   motion_cache : containers.Map, fly id -> per-frame sum_motion vector.
%
%   opts (struct, all optional — defaults in local geto()):
%       fps(60) kernel_past_s(3) kernel_future_s(1.5) dt_frames(6)
%       nb_kernel(6) nb_acausal(6) bump_width(1.5) nb_baseline(6)
%       cv_folds(5) n_shuffle(0) shuffle_mode('full') trunc_point(30)
%       do_plot(false) do_stage2(true) use_penalty(true) fig_title('') rng_seed(7)
%
%   use_penalty: P-spline 2nd-difference roughness penalty on the kernel coeffs,
%       λ chosen by the 1-SE rule (damps ringing from collinear/autocorrelated
%       lags). false → plain unpenalised logistic kernel.
%
%   res fields: kernel, kernel_se, kernel_unpen, t_rel, lag_fr, caus_mask,
%       beta_k, Cov_k, Bk, lambda, lam_grid, cvll_lam, cvse_lam,
%       tau_hat, tau_ci, R2_exp, tau_com, aic, cvll, specs, n_events, n_int,
%       acf_sm, eta_ev, eta_exc, rho_acf, rho_eta (kernel-shape diagnostics),
%       and (if n_shuffle>0) ll_real, ll_null, pval, kernels_sh.

% ── options ──────────────────────────────────────────────────────────────
if nargin < 3, opts = struct(); end
fps        = geto(opts,'fps',60);
kpast      = geto(opts,'kernel_past_s',3);
kfut       = geto(opts,'kernel_future_s',3);
dt_frames  = geto(opts,'dt_frames',6);
mask_preonset = geto(opts,'mask_preonset',false);  % zero out pre-onset lags
grid_anchor = geto(opts,'grid_anchor','entry');    % 'entry' from onset | 'offset' from break
nb_kernel  = geto(opts,'nb_kernel',6);
nb_acausal = geto(opts,'nb_acausal',6);
bump_width = geto(opts,'bump_width',1.5);
nb_baseline= geto(opts,'nb_baseline',6);
cv_folds   = geto(opts,'cv_folds',5);
n_shuffle  = geto(opts,'n_shuffle',0);
shuffle_mode = geto(opts,'shuffle_mode','full');
trunc_point  = geto(opts,'trunc_point',30);
do_plot    = geto(opts,'do_plot',false);
do_stage2  = geto(opts,'do_stage2',true);
do_unpen   = geto(opts,'do_unpen',true);       % compute the unpenalised overlay kernel (plain
                                               % fitglm). Skip for big/collinear designs where the
                                               % MLE may not converge — only the PENALISED kernel
                                               % is reported, so this is purely cosmetic.
use_penalty= geto(opts,'use_penalty',true);   % P-spline roughness penalty on the kernel (λ by 1-SE)
kernel_ridge = geto(opts,'kernel_ridge',1e-3); % small λ-independent ridge on the kernel coeffs (×kernel
                                               % Gram). Removes the roughness null space so the kernel
                                               % amplitude can't diverge when the likelihood is flat
                                               % (e.g. circshift-null fits with the sm signal removed).
lam_fixed  = geto(opts,'lambda_fixed',[]);     % fix λ (skip the CV grid); used by the circular-shift
                                               % null so every shuffle shares the real-fit smoothing
do_leakfit = geto(opts,'do_leakfit',true);     % leak-τ exponential fit + 2000× bootstrap CI; set false
                                               % when only the kernel is needed (e.g. the circshift null)
circ_min_shift_s = geto(opts,'circ_min_shift_s',5);   % 'circshift' mode: min |rotation| per sm trace (s)
anchor_zero= geto(opts,'anchor_zero',true);    % 'full' mode: force a bump center on τ=0
fig_title  = geto(opts,'fig_title','');
rng_seed   = geto(opts,'rng_seed',7);
rng(rng_seed);

entry_fr = trunc_point;

% ── lag grid + kernel basis ───────────────────────────────────────────────
back_fr =  round(kpast*fps);
fwd_fr  =  round(kfut *fps);
lag_fr  = (-fwd_fr:back_fr)';
nLag    = numel(lag_fr);
caus_mask = lag_fr >= 0;
caus_t    = lag_fr(caus_mask)/fps;

anchor_lag = []; if anchor_zero, anchor_lag = 0; end   % knot on τ=0 (lag 0) when requested
switch shuffle_mode
    case {'full','circshift'}
        Bk = raised_cosine_basis(lag_fr, nb_kernel+nb_acausal, bump_width, anchor_lag);
    case 'split'
        Bk = zeros(nLag, nb_kernel+nb_acausal);
        Bk(caus_mask,  1:nb_kernel)              = raised_cosine_basis(lag_fr(caus_mask),  nb_kernel,  bump_width);
        Bk(~caus_mask, nb_kernel+(1:nb_acausal)) = raised_cosine_basis(lag_fr(~caus_mask), nb_acausal, bump_width);
    otherwise
        error('shuffle_mode must be ''full'', ''split'', or ''circshift''');
end
nb_k = size(Bk,2);

% ── global z-scoring of sm (+ cache traces for the circshift null) ─────────
is_circ    = strcmp(shuffle_mode,'circshift');
keep_shift = n_shuffle > 0 && is_circ;          % store per-row fly/frame + full traces
flies = unique(bl.fly);  nF = numel(flies);
[~, bl_flyidx] = ismember(bl.fly, flies);       % bout -> fly index (1..nF)
smcell = cell(nF,1); Lf = zeros(nF,1);
nn=0; s1=0; s2=0;
for i=1:nF
    v = motion_cache(flies(i)); v = v(:);
    if keep_shift, smcell{i} = v; Lf(i) = numel(v); end
    vv = v(~isnan(v));
    nn=nn+numel(vv); s1=s1+sum(vv); s2=s2+sum(vv.^2);
end
sm_mu = s1/nn;  sm_sd = sqrt(s2/nn - sm_mu^2);

% ── build long (person-period) design ─────────────────────────────────────
max_rows = sum(ceil(bl.dur_frames/dt_frames)) + height(bl);
Xk=nan(max_rows,nb_k); tinf=nan(max_rows,1);
yev=zeros(max_rows,1); gid=nan(max_rows,1); mov=nan(max_rows,1); slo=nan(max_rows,1);
fmean=nan(max_rows,1); fmax=nan(max_rows,1); fslp=nan(max_rows,1);

% The full per-interval sm history S (max_rows x nLag) is ONLY needed for the
% time-shuffle control; for large/augmented datasets it can be tens of GB, so
% allocate it only when shuffling. The ACF / event-triggered diagnostics are
% accumulated online below, so they never need the full matrix.
keep_S = n_shuffle > 0 && ~is_circ;          % full/split need S; circshift re-reads shifted traces
if keep_S, S = nan(max_rows, nLag); else, S = []; end
if keep_shift, tabs = nan(max_rows,1); frow = nan(max_rows,1); end
zero_idx = find(lag_fr==0, 1);
acc_s = zeros(nLag,1); acc_ss0 = zeros(nLag,1); acc_se = zeros(nLag,1); n_ev = 0;

A_slope = caus_t - mean(caus_t);
A_slope = A_slope / (A_slope'*A_slope);

r = 0;
for b = 1:height(bl)
    sm = motion_cache(bl.fly(b)); sm = sm(:); L = numel(sm);
    on = bl.onsets(b); dur = bl.dur_frames(b);
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
        if (t_abs-back_fr) < 1 || (t_abs+fwd_fr) > L, continue; end
        s = (sm(t_abs - lag_fr) - sm_mu)/sm_sd;
        if any(isnan(s)), continue; end
        if mask_preonset, s(lag_fr > f-1) = 0; end   % drop pre-onset drive
        r = r+1;
        Xk(r,:) = s'*Bk;
        if keep_S, S(r,:) = s'; end
        if keep_shift, tabs(r) = t_abs; frow(r) = bl_flyidx(b); end
        tinf(r) = f/fps;   gid(r) = bl.bout_id(b);
        mov(r)  = bl.moving_flies(b);  slo(r) = bl.sloom(b);
        sc = s(caus_mask);
        fmean(r)=mean(sc); fmax(r)=max(sc); fslp(r)=A_slope'*(sc-mean(sc));
        % online accumulators for the ACF / event-triggered diagnostics
        acc_s = acc_s + s;  acc_ss0 = acc_ss0 + s*s(zero_idx);
        if gi==numel(grid) && ~bl.is_censored(b)
            yev(r)=1; acc_se = acc_se + s; n_ev = n_ev + 1;
        end
    end
end
Xk=Xk(1:r,:); tinf=tinf(1:r); yev=yev(1:r); gid=gid(1:r);
if keep_S, S=S(1:r,:); end
if keep_shift, tabs=tabs(1:r); frow=frow(1:r); end
mov=mov(1:r); slo=slo(1:r); fmean=fmean(1:r); fmax=fmax(1:r); fslp=fslp(1:r);
res.n_int = r; res.n_events = sum(yev);

% ── baseline-hazard basis (rank-scale knots, drop one for intercept) ──────
ut = tiedrank(tinf)/numel(tinf);
Bb = raised_cosine_basis(ut, nb_baseline);
Bb = Bb(:,2:end);
nb_base_use = size(Bb,2);

% ── controls: include mov/slo only if they vary ──────────────────────────
ctrl_tbl = table();  ctrl_terms = '';  ctrl_num = [];
if numel(unique(mov))>1
    ctrl_tbl.mov = categorical(mov); ctrl_terms=[ctrl_terms ' + mov'];
    [~,~,im]=unique(mov);  Dm=dummyvar(im);  ctrl_num=[ctrl_num, Dm(:,2:end)];
end
if numel(unique(slo))>1
    ctrl_tbl.slo = categorical(slo); ctrl_terms=[ctrl_terms ' + slo'];
    [~,~,iss]=unique(slo); Ds=dummyvar(iss); ctrl_num=[ctrl_num, Ds(:,2:end)];
end
hterms = strjoin(arrayfun(@(j)sprintf('h%d',j),1:nb_base_use,'uni',0),' + ');

% fold assignment (grouped by bout) — needed by λ-selection, Stage 2, shuffle
fold = mod(randperm(max(gid)),cv_folds)+1; rowfold = fold(gid);

% ── numeric design + P-spline roughness penalty on the kernel coeffs ──────
% Mirrors hazard_kernel_sm.m: a 2nd-difference penalty on the evenly-spaced
% raised-cosine weights damps the +,−,+ ringing from collinear, autocorrelated
% lags; λ by the 1-SE rule. use_penalty=false → plain unpenalised kernel.
Xdes = [ones(r,1), Xk, Bb, ctrl_num];  ker_cols = 1 + (1:nb_k);
switch shuffle_mode
    case {'full','circshift'}, Dk = diffmat(nb_k,2);
    case 'split',              Dk = blkdiag(diffmat(nb_kernel,2), diffmat(nb_acausal,2));
end
Pk = zeros(size(Xdes,2)); Pk(ker_cols,ker_cols) = Dk'*Dk;
% λ-independent ridge on the kernel coeffs. Dk'Dk has a null space ({flat, linear}
% kernels), so when the likelihood is flat along the kernel directions (e.g. a
% circshift-null fit where sm no longer predicts escapes) the amplitude is
% unidentified and can diverge. A small ridge — scaled to the kernel Gram so it is
% data-adaptive, and NOT multiplied by λ so it survives a small λ — removes the
% null space and bounds the amplitude; it is negligible vs the likelihood curvature
% when real signal is present. Added to BOTH the CV objective and the final fit.
Pr = zeros(size(Xdes,2));
ridge_abs = kernel_ridge * mean(diag(Xk'*Xk));
Pr(ker_cols,ker_cols) = ridge_abs * eye(nb_k);
res.ridge_abs = ridge_abs;

if use_penalty
    if ~isempty(lam_fixed)
        % λ pinned by the caller (e.g. the circular-shift null) — skip the CV grid.
        lambda = lam_fixed; lam_grid = lam_fixed; cvll_lam = NaN; cvse_lam = NaN;
    else
        lam_grid = logspace(-1,6,15); cvll_lam=nan(size(lam_grid)); cvse_lam=nan(size(lam_grid));
        for il=1:numel(lam_grid)
            [cvll_lam(il), llf] = pen_grouped_cv(Xdes, yev, lam_grid(il)*Pk + Pr, rowfold, cv_folds);
            cvse_lam(il) = std(llf,'omitnan')/sqrt(sum(~isnan(llf)));
        end
        [~,ibest] = max(cvll_lam);
        i1se   = find(cvll_lam >= cvll_lam(ibest)-cvse_lam(ibest), 1, 'last');  % 1-SE rule
        lambda = lam_grid(i1se);
    end
    P_pen = lambda*Pk + Pr;
else
    lam_grid=0; cvll_lam=NaN; cvse_lam=NaN; lambda=0;
    P_pen = lambda*Pk;             % truly unpenalised (no ridge)
end

% ── kernel: unpenalised (fitglm, overlay) + penalised (reported) ──────────
% The unpenalised MLE is only a cosmetic overlay; on large/collinear designs it
% may hit glmfit's iteration cap (a benign warning). Skip it when do_unpen=false.
if do_unpen
    [mdl, ~] = fit_kernel_glm(Xk, Bb, ctrl_tbl, ctrl_terms, yev);
    cn = mdl.CoefficientNames; kidx = zeros(nb_k,1);
    for j=1:nb_k, kidx(j)=find(strcmp(cn,sprintf('k%d',j))); end
    kernel_unpen = Bk*mdl.Coefficients.Estimate(kidx);
else
    kernel_unpen = nan(size(Bk,1),1);
end

[beta_all, V_all] = penalized_logit(Xdes, yev, P_pen);
beta_k = beta_all(ker_cols);
Cov_k  = V_all(ker_cols,ker_cols);
kernel    = Bk*beta_k;
kernel_se = sqrt(sum((Bk*Cov_k).*Bk,2));
t_rel     = -lag_fr/fps;

res.kernel=kernel; res.kernel_se=kernel_se; res.kernel_unpen=kernel_unpen;
res.t_rel=t_rel; res.lag_fr=lag_fr; res.caus_mask=caus_mask;
res.beta_k=beta_k; res.Cov_k=Cov_k; res.Bk=Bk;
res.use_penalty=use_penalty; res.lambda=lambda;
res.lam_grid=lam_grid; res.cvll_lam=cvll_lam; res.cvse_lam=cvse_lam;

% ── kernel-shape diagnostic: sm autocorrelation + event-triggered average ──
% A symmetric β(τ) is the signature of an autocorrelated regressor (sm elevated
% AROUND escape) rather than directed causality. Store R(τ), the event-triggered
% sm excess, and how strongly the kernel correlates with each.
% Computed from the online accumulators (no full-S storage needed).
if n_ev > 0, eta_ev = acc_se / n_ev; else, eta_ev = nan(nLag,1); end
mean_s  = acc_s / max(r,1);
eta_exc = eta_ev - mean_s;
acf_sm  = acc_ss0 / max(r,1);  acf_sm = acf_sm / acf_sm(zero_idx);
res.acf_sm=acf_sm; res.eta_ev=eta_ev; res.eta_exc=eta_exc;
res.rho_acf = corr(kernel, acf_sm);
res.rho_eta = corr(kernel, eta_exc, 'rows','complete');

% ── leak time-constant: exponential fit + model-free COM + bootstrap CI ───
% Skipped when do_leakfit=false (e.g. the circshift null only needs the kernel) —
% the 2000× parametric bootstrap is the single most expensive post-fit step.
if do_leakfit
    cm=caus_mask; tau_s=lag_fr(cm)/fps; kc=kernel(cm); sec=kernel_se(cm);
    wts = 1./max(sec,eps).^2;
    expmod = @(p,t) p(1).*exp(-t./exp(p(2))) + p(3);
    sse = @(p) sum(wts.*(kc-expmod(p,tau_s)).^2);
    p0 = [max(kc)-min(kc), log(0.5), min(kc)];
    oopt = optimset('MaxFunEvals',1e4,'MaxIter',1e4,'Display','off');
    pfit = fminsearch(sse, p0, oopt);
    tau_hat = exp(pfit(2));
    resid = kc-expmod(pfit,tau_s); wm = sum(wts.*kc)/sum(wts);
    R2_exp = 1 - sum(wts.*resid.^2)/sum(wts.*(kc-wm).^2);
    kpos = max(kc,0);
    if sum(kpos)>0, tau_com = sum(tau_s.*kpos)/sum(kpos); else, tau_com = NaN; end
    % parametric bootstrap CI
    nboot=2000; Csym=(Cov_k+Cov_k')/2; [Vc,Dc]=eig(Csym); Csym=Vc*diag(max(diag(Dc),0))*Vc';
    Bsamp = mvnrnd(beta_k(:)',Csym,nboot); Bk_c=Bk(cm,:); tb_all=nan(nboot,1);
    for ii=1:nboot
        kb=Bk_c*Bsamp(ii,:)';
        try
            pb=fminsearch(@(p)sum(wts.*(kb-expmod(p,tau_s)).^2),pfit,oopt);
            tt=exp(pb(2));
            if tt>1e-3 && tt<max(tau_s)*5, tb_all(ii)=tt; end
        catch
        end
    end
    tau_ci = prctile(tb_all(~isnan(tb_all)),[2.5 97.5]);
    res.tau_hat=tau_hat; res.tau_ci=tau_ci; res.R2_exp=R2_exp; res.tau_com=tau_com;
else
    res.tau_hat=NaN; res.tau_ci=[NaN NaN]; res.R2_exp=NaN; res.tau_com=NaN;
end

% (fold assignment `rowfold` already done above, before λ-selection)

% ── Stage 2: integration vs extrema model comparison (optional; slow) ─────
res.aic=[]; res.cvll=[]; res.specs=[];
if do_stage2
    zc=@(x)(x-mean(x))./std(x);
    S2 = ctrl_tbl; S2.y = yev;
    for j=1:nb_base_use, S2.(sprintf('h%d',j))=Bb(:,j); end
    S2.f_mean=zc(fmean); S2.f_max=zc(fmax); S2.f_slope=zc(fslp);
    ctrl = ['1 + ' hterms ctrl_terms];
    specs = struct('name',{'baseline','mean','max','slope','mean+max','mean+slope','mean+max+slope'}, ...
        'rhs',{ctrl, ['f_mean + ' ctrl], ['f_max + ' ctrl], ['f_slope + ' ctrl], ...
               ['f_mean + f_max + ' ctrl], ['f_mean + f_slope + ' ctrl], ['f_mean + f_max + f_slope + ' ctrl]});
    aic=nan(numel(specs),1); cvll=nan(numel(specs),1);
    for m=1:numel(specs)
        mm=fitglm(S2,['y ~ ' specs(m).rhs],'Distribution','binomial','Options',statset('MaxIter',1000));
        aic(m)=mm.ModelCriterion.AIC;
        cvll(m)=grouped_cv_loglik(S2,['y ~ ' specs(m).rhs],rowfold,cv_folds);
    end
    res.aic=aic; res.cvll=cvll; res.specs=specs;
end

% ── optional shuffle control ──────────────────────────────────────────────
if n_shuffle > 0
    Xdes_sh = Xdes;                              % template; only kernel cols change
    kernels_sh = nan(nLag, n_shuffle);
    if is_circ
        % ── circular-shift null ────────────────────────────────────────────
        % Rotate every fly's sm trace by an independent random offset: preserves
        % its autocorrelation, breaks alignment to escapes. The design is built
        % ONCE; each shuffle only re-reads the shifted windows (vectorised per
        % fly) and refits the penalised logit warm-started from the real β at the
        % fixed λ — no rebuild, no CV, no bootstrap. Any window straddling the
        % wrap seam has its kernel regressor blanked to the mean (0 after z-score).
        span = back_fr + fwd_fr;
        min_shift = max(round(circ_min_shift_s*fps), span+1);
        rows_fi = cell(nF,1); tabs_fi = cell(nF,1);
        for fi=1:nF, rr=find(frow==fi); rows_fi{fi}=rr; tabs_fi{fi}=tabs(rr); end
        nseam = 0;  Ssh = zeros(r, nLag);        % preallocated; every row overwritten per shuffle
        for sN=1:n_shuffle
            for fi=1:nF
                Lfi = Lf(fi); if isempty(rows_fi{fi}), continue; end
                if Lfi > 2*min_shift, d = randi([min_shift, Lfi-min_shift]); else, d = max(1,floor(Lfi/2)); end
                raw = tabs_fi{fi} + d - lag_fr.';             % (n_fi x nLag) absolute positions
                wrapped = floor((min(raw,[],2)-1)/Lfi) ~= floor((max(raw,[],2)-1)/Lfi);
                vals = (smcell{fi}(mod(raw-1,Lfi)+1) - sm_mu)/sm_sd;
                vals(isnan(vals)) = 0;  vals(wrapped,:) = 0;  % NaN gaps + seam -> mean
                Ssh(rows_fi{fi},:) = vals;
                if sN==1, nseam = nseam + sum(wrapped); end
            end
            Xdes_sh(:,ker_cols) = Ssh*Bk;
            b_sh = penalized_logit(Xdes_sh, yev, P_pen, beta_all);   % warm start from real β
            kernels_sh(:,sN) = Bk*b_sh(ker_cols);
        end
        res.kernels_sh=kernels_sh; res.ll_real=NaN; res.ll_null=nan(n_shuffle,1);
        res.pval=NaN; res.seam_frac=nseam/max(r,1);
    else
        % ── within-window temporal-order shuffle ('full' / 'split') ─────────
        caus_idx=find(caus_mask); nca=numel(caus_idx);
        acaus_idx=find(~caus_mask); naa=numel(acaus_idx);
        ll_real = pen_grouped_cv(Xdes, yev, P_pen, rowfold, cv_folds);
        ll_null=nan(n_shuffle,1);
        nr=size(S,1); rows_all=repmat((1:nr)',1,nLag);
        Ac=S(:,caus_idx); rows_c=repmat((1:nr)',1,nca);
        Aa=S(:,acaus_idx); rows_a=repmat((1:nr)',1,naa);
        for sN=1:n_shuffle
            switch shuffle_mode
                case 'full'
                    [~,ord]=sort(rand(nr,nLag),2);
                    Ssh=S(sub2ind([nr nLag],rows_all,ord));
                case 'split'
                    [~,ordc]=sort(rand(nr,nca),2); [~,orda]=sort(rand(nr,naa),2);
                    Ssh=S; Ssh(:,caus_idx)=Ac(sub2ind([nr nca],rows_c,ordc));
                    Ssh(:,acaus_idx)=Aa(sub2ind([nr naa],rows_a,orda));
            end
            Xdes_sh(:,ker_cols)=Ssh*Bk;
            ll_null(sN)=pen_grouped_cv(Xdes_sh,yev,P_pen,rowfold,cv_folds);
            b_sh=penalized_logit(Xdes_sh,yev,P_pen);
            kernels_sh(:,sN)=Bk*b_sh(ker_cols);
        end
        res.ll_real=ll_real; res.ll_null=ll_null;
        res.pval=(1+sum(ll_null>=ll_real))/(1+n_shuffle); res.kernels_sh=kernels_sh;
    end
end

% ── optional plot ─────────────────────────────────────────────────────────
if do_plot
    figure('Color','w','Position',[100 100 720 460]); hold on
    fill([t_rel;flipud(t_rel)],[kernel+1.96*kernel_se;flipud(kernel-1.96*kernel_se)], ...
        [0.2 0.5 0.8],'FaceAlpha',0.2,'EdgeColor','none')
    plot(t_rel,kernel,'Color',[0.2 0.5 0.8],'LineWidth',2)
    yline(0,'k:'); xline(0,'--k','escape')
    xlabel('Time relative to un-freeze (s)   —   past < 0 < future')
    ylabel('Kernel weight  \beta(\tau)')
    title(sprintf('%s   (\\tau_{com}=%.2f s, \\tau_{exp}=%.2f s, R²=%.2f)', ...
        fig_title, tau_com, tau_hat, R2_exp))
end
end

% ════════════════════════════════════════════════════════════════════════
function v = geto(s,f,d), if isfield(s,f) && ~isempty(s.(f)), v=s.(f); else, v=d; end, end

function [mdl, T] = fit_kernel_glm(Xk_in, Bb, ctrl_tbl, ctrl_terms, yev)
    T = ctrl_tbl; T.y = yev;
    nk=size(Xk_in,2); nb=size(Bb,2);
    for j=1:nk, T.(sprintf('k%d',j))=Xk_in(:,j); end
    for j=1:nb, T.(sprintf('h%d',j))=Bb(:,j); end
    kt=strjoin(arrayfun(@(j)sprintf('k%d',j),1:nk,'uni',0),' + ');
    ht=strjoin(arrayfun(@(j)sprintf('h%d',j),1:nb,'uni',0),' + ');
    mdl=fitglm(T, sprintf('y ~ %s + %s%s', kt, ht, ctrl_terms), 'Distribution','binomial','Link','logit', ...
        'Options', statset('MaxIter', 1000));   % default glmfit cap is 100; raise it so the (collinear) overlay converges
end

function D = diffmat(n, ord)
% k-th order difference operator: (n-ord) x n. diffmat(n,2) penalises curvature.
    D = eye(n);
    for k = 1:ord, D = diff(D); end
end

function [b, V] = penalized_logit(X, y, P, b0)
% Penalised logistic regression by IRLS/Newton: maximise loglik − 0.5·bᵀPb.
% P is the (p×p) quadratic penalty (zeros for unpenalised coeffs). V = inv(XᵀWX+P)
% is the P-spline posterior covariance for the kernel CI.
% b0 (optional): warm-start coefficients (e.g. the real-fit β for shuffle refits).
    [~, p] = size(X);  rdg = 1e-8*eye(p);
    if nargin >= 4 && ~isempty(b0), b = b0(:); else, b = zeros(p,1); end
    for it = 1:100
        eta = min(max(X*b,-30),30);
        mu  = 1./(1+exp(-eta));
        w   = max(mu.*(1-mu), 1e-9);
        H   = X'*(X.*w) + P + rdg;
        step = H \ (X'*(y-mu) - P*b);
        if ~all(isfinite(step)), break; end          % singular/divergent step -> stop
        sl = max(abs(step));                          % damped Newton: cap the per-iter step so a
        if sl > 10, step = step*(10/sl); end          % near-singular H can't launch b to ±∞
        b = b + step;
        if sl < 1e-8, break; end
    end
    V = inv(H);
end

function [ll, ll_folds] = pen_grouped_cv(X, y, P, rowfold, K)
% Held-out Bernoulli logLik/obs of the penalised model, folds grouped by bout.
% Returns across-fold mean and per-fold means (the latter feeds the 1-SE rule).
    ll_folds = nan(K,1);
    for k = 1:K
        te = rowfold==k; tr = ~te;
        if ~any(te)||~any(tr), continue; end
        b   = penalized_logit(X(tr,:), y(tr), P);
        eta = min(max(X(te,:)*b,-30),30);
        mu  = min(max(1./(1+exp(-eta)),1e-9),1-1e-9);
        yy  = y(te);
        ll_folds(k) = mean(yy.*log(mu) + (1-yy).*log(1-mu));
    end
    ll = mean(ll_folds,'omitnan');
end

function ll = grouped_cv_loglik(T, formula, rowfold, K)
    acc=0; cnt=0;
    for k=1:K
        te=rowfold==k; tr=~te;
        if ~any(te)||~any(tr), continue; end
        mm=fitglm(T(tr,:),formula,'Distribution','binomial','Options',statset('MaxIter',1000));
        p=predict(mm,T(te,:)); p=min(max(p,1e-9),1-1e-9); y=T.y(te);
        acc=acc+sum(y.*log(p)+(1-y).*log(1-p)); cnt=cnt+sum(te);
    end
    ll=acc/cnt;
end

function B = raised_cosine_basis(x, nb, width_mult, anchor)
% anchor (optional): shift the evenly-spaced grid by <½ spacing so one center lands
%   exactly on `anchor` (e.g. 0 → a knot at τ=0). Spacing stays uniform.
    if nargin<3||isempty(width_mult), width_mult=1; end
    if nargin<4, anchor=[]; end
    x=x(:); lo=min(x); hi=max(x);
    c=linspace(lo,hi,nb); sp=(hi-lo)/(nb-1); w=width_mult*sp;
    if ~isempty(anchor)
        [~,jn]=min(abs(c-anchor)); c=c+(anchor-c(jn));
    end
    B=zeros(numel(x),nb);
    for j=1:nb
        d=(x-c(j))/w; in=abs(d)<=1;
        B(in,j)=0.5*(1+cos(pi*d(in)));
    end
end
