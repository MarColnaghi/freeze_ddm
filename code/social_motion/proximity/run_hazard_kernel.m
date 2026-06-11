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
use_penalty= geto(opts,'use_penalty',true);   % P-spline roughness penalty on the kernel (λ by 1-SE)
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
    case 'full'
        Bk = raised_cosine_basis(lag_fr, nb_kernel+nb_acausal, bump_width, anchor_lag);
    case 'split'
        Bk = zeros(nLag, nb_kernel+nb_acausal);
        Bk(caus_mask,  1:nb_kernel)              = raised_cosine_basis(lag_fr(caus_mask),  nb_kernel,  bump_width);
        Bk(~caus_mask, nb_kernel+(1:nb_acausal)) = raised_cosine_basis(lag_fr(~caus_mask), nb_acausal, bump_width);
    otherwise
        error('shuffle_mode must be ''full'' or ''split''');
end
nb_k = size(Bk,2);

% ── global z-scoring of sm ────────────────────────────────────────────────
flies = unique(bl.fly);
nn=0; s1=0; s2=0;
for i=1:numel(flies)
    v = motion_cache(flies(i)); v = v(:); v = v(~isnan(v));
    nn=nn+numel(v); s1=s1+sum(v); s2=s2+sum(v.^2);
end
sm_mu = s1/nn;  sm_sd = sqrt(s2/nn - sm_mu^2);

% ── build long (person-period) design ─────────────────────────────────────
max_rows = sum(ceil(bl.dur_frames/dt_frames)) + height(bl);
Xk=nan(max_rows,nb_k); tinf=nan(max_rows,1);
S=nan(max_rows,nLag);   % raw z-scored history (shuffle control + ACF/ETA diagnostic)
yev=zeros(max_rows,1); gid=nan(max_rows,1); mov=nan(max_rows,1); slo=nan(max_rows,1);
fmean=nan(max_rows,1); fmax=nan(max_rows,1); fslp=nan(max_rows,1);

A_slope = caus_t - mean(caus_t);
A_slope = A_slope / (A_slope'*A_slope);

r = 0;
for b = 1:height(bl)
    sm = motion_cache(bl.fly(b)); sm = sm(:); L = numel(sm);
    on = bl.onsets(b); dur = bl.dur_frames(b);
    if dur < entry_fr, continue; end
    grid = (entry_fr:dt_frames:dur)';
    if isempty(grid) || grid(end) < dur, grid = [grid; dur]; end %#ok<AGROW>
    for gi = 1:numel(grid)
        f = grid(gi); t_abs = on + f - 1;
        if (t_abs-back_fr) < 1 || (t_abs+fwd_fr) > L, continue; end
        s = (sm(t_abs - lag_fr) - sm_mu)/sm_sd;
        if any(isnan(s)), continue; end
        r = r+1;
        Xk(r,:) = s'*Bk;   S(r,:) = s';
        tinf(r) = f/fps;   gid(r) = bl.bout_id(b);
        mov(r)  = bl.moving_flies(b);  slo(r) = bl.sloom(b);
        sc = s(caus_mask);
        fmean(r)=mean(sc); fmax(r)=max(sc); fslp(r)=A_slope'*(sc-mean(sc));
        if gi==numel(grid) && ~bl.is_censored(b), yev(r)=1; end
    end
end
Xk=Xk(1:r,:); tinf=tinf(1:r); yev=yev(1:r); gid=gid(1:r); S=S(1:r,:);
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
    case 'full',  Dk = diffmat(nb_k,2);
    case 'split', Dk = blkdiag(diffmat(nb_kernel,2), diffmat(nb_acausal,2));
end
Pk = zeros(size(Xdes,2)); Pk(ker_cols,ker_cols) = Dk'*Dk;

if use_penalty
    lam_grid = logspace(-1,6,15); cvll_lam=nan(size(lam_grid)); cvse_lam=nan(size(lam_grid));
    for il=1:numel(lam_grid)
        [cvll_lam(il), llf] = pen_grouped_cv(Xdes, yev, lam_grid(il)*Pk, rowfold, cv_folds);
        cvse_lam(il) = std(llf,'omitnan')/sqrt(sum(~isnan(llf)));
    end
    [~,ibest] = max(cvll_lam);
    i1se   = find(cvll_lam >= cvll_lam(ibest)-cvse_lam(ibest), 1, 'last');  % 1-SE rule
    lambda = lam_grid(i1se);
else
    lam_grid=0; cvll_lam=NaN; cvse_lam=NaN; lambda=0;
end
P_pen = lambda*Pk;

% ── kernel: unpenalised (fitglm, overlay) + penalised (reported) ──────────
[mdl, ~] = fit_kernel_glm(Xk, Bb, ctrl_tbl, ctrl_terms, yev);
cn = mdl.CoefficientNames; kidx = zeros(nb_k,1);
for j=1:nb_k, kidx(j)=find(strcmp(cn,sprintf('k%d',j))); end
kernel_unpen = Bk*mdl.Coefficients.Estimate(kidx);

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
zero_idx = find(lag_fr==0,1); ev = yev==1;
if any(ev), eta_ev = mean(S(ev,:),1)'; else, eta_ev = nan(nLag,1); end
eta_exc = eta_ev - mean(S,1)';
acf_sm  = mean(S.*S(:,zero_idx),1)'; acf_sm = acf_sm/acf_sm(zero_idx);
res.acf_sm=acf_sm; res.eta_ev=eta_ev; res.eta_exc=eta_exc;
res.rho_acf = corr(kernel, acf_sm);
res.rho_eta = corr(kernel, eta_exc, 'rows','complete');

% ── leak time-constant: exponential fit + model-free COM + bootstrap CI ───
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
        mm=fitglm(S2,['y ~ ' specs(m).rhs],'Distribution','binomial');
        aic(m)=mm.ModelCriterion.AIC;
        cvll(m)=grouped_cv_loglik(S2,['y ~ ' specs(m).rhs],rowfold,cv_folds);
    end
    res.aic=aic; res.cvll=cvll; res.specs=specs;
end

% ── optional time-shuffle control ─────────────────────────────────────────
if n_shuffle > 0
    caus_idx=find(caus_mask); nca=numel(caus_idx);
    acaus_idx=find(~caus_mask); naa=numel(acaus_idx);
    ll_real = pen_grouped_cv(Xdes, yev, P_pen, rowfold, cv_folds);
    ll_null=nan(n_shuffle,1); kernels_sh=nan(nLag,n_shuffle);
    nr=size(S,1); rows_all=repmat((1:nr)',1,nLag);
    Ac=S(:,caus_idx); rows_c=repmat((1:nr)',1,nca);
    Aa=S(:,acaus_idx); rows_a=repmat((1:nr)',1,naa);
    Xdes_sh = Xdes;                              % template; only kernel cols re-shuffled
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
    mdl=fitglm(T, sprintf('y ~ %s + %s%s', kt, ht, ctrl_terms), 'Distribution','binomial','Link','logit');
end

function D = diffmat(n, ord)
% k-th order difference operator: (n-ord) x n. diffmat(n,2) penalises curvature.
    D = eye(n);
    for k = 1:ord, D = diff(D); end
end

function [b, V] = penalized_logit(X, y, P)
% Penalised logistic regression by IRLS/Newton: maximise loglik − 0.5·bᵀPb.
% P is the (p×p) quadratic penalty (zeros for unpenalised coeffs). V = inv(XᵀWX+P)
% is the P-spline posterior covariance for the kernel CI.
    [~, p] = size(X);  b = zeros(p,1);  rdg = 1e-8*eye(p);
    for it = 1:100 %#ok<NASGU>
        eta = min(max(X*b,-30),30);
        mu  = 1./(1+exp(-eta));
        w   = max(mu.*(1-mu), 1e-9);
        H   = X'*(X.*w) + P + rdg;
        step = H \ (X'*(y-mu) - P*b);
        b = b + step;
        if max(abs(step)) < 1e-8, break; end
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
        mm=fitglm(T(tr,:),formula,'Distribution','binomial');
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
