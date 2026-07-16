function out = fit_kernels_sparse(S, yev, tinf, mov, slo, gid, lag_fr, fps, opts)
% FIT_KERNELS_SPARSE  Estimate the un-freeze hazard kernel β(τ) from a person-period
% design, four ways, following Mineault, Barthelmé & Pack (2009). This is the
% function-form core of hazard_kernel_sparse.m so that the real analysis and the
% recovery experiments run IDENTICAL estimation code.
%
%   S      : (r × nLag) raw z-scored sm history per at-risk interval
%   yev    : (r × 1) 1 = un-freeze event this interval
%   tinf   : (r × 1) time-in-freeze (s)  -> baseline hazard
%   mov,slo: (r × 1) covariates (controls)
%   gid    : (r × 1) bout id (grouped CV)
%   lag_fr : (nLag × 1) lag grid in frames (τ = lag/fps; ≥0 causal, <0 acausal)
%   fps    : frames/s
%   opts   : struct (all optional): nb_kernel, nb_acausal, bump_width, nb_baseline,
%            n_exp, cv_folds, kernel_past_s, sel_rule ('argmax'|'1se'), verbose
%
% Returns struct OUT with: t_rel, caus_mask, lag_fr, fps, tau_grid, is_exp, sel_tau,
%   kernel_A/B/C (+_se), per-model dev/df/aic/cv_dev/lambda, and the nuisance-only
%   baseline (suffix N). Models: (0) nuisance-only, (A) L2 smoothness on raised
%   cosines, (B) L1 sparse on an exponential bank, (C) L1 sparse on flexible bumps.

    o = default_opts(opts);
    r = size(S,1);
    caus_mask = lag_fr >= 0;
    t_rel     = -lag_fr / fps;

    % ── baseline-hazard basis (equal-mass knots; drop one) + nuisance dummies ──
    ut = tiedrank(tinf) / numel(tinf);
    Bb = raised_cosine_basis(ut, o.nb_baseline);  Bb = Bb(:, 2:end);
    [~,~,im]  = unique(mov);  Dmov = dummyvar(im);  Dmov = Dmov(:,2:end);
    [~,~,iss] = unique(slo);  Dslo = dummyvar(iss); Dslo = Dslo(:,2:end);
    Nuis = [Bb, Dmov, Dslo];

    fold = mod(randperm(max(gid)), o.cv_folds) + 1;  rowfold = fold(gid);

    % ── the three kernel bases ────────────────────────────────────────────────
    tau_grid_s = logspace(log10(0.05), log10(o.kernel_past_s), o.n_exp);
    [Bexp, tau_grid, is_exp] = exp_basis(lag_fr, tau_grid_s, fps, o.nb_acausal, o.bump_width);
    Bflex = raised_cosine_basis(lag_fr, o.nb_kernel + o.nb_acausal, o.bump_width);
    Bflex = Bflex ./ max(vecnorm(Bflex), eps);
    Bsm   = raised_cosine_basis(lag_fr, o.nb_kernel + o.nb_acausal, o.bump_width);

    % ── (0) nuisance-only baseline ────────────────────────────────────────────
    Xnull = [ones(r,1), Nuis];  PtinyN = 1e-6 * eye(size(Xnull,2));
    [bN, VN] = penalized_logit(Xnull, yev, PtinyN);
    [devN, dfN, aicN] = model_metrics_l2(Xnull, yev, bN, VN);
    cv_devN = -2 * pen_grouped_cv(Xnull, yev, PtinyN, rowfold, o.cv_folds);

    % ── (A) L2 smoothness ─────────────────────────────────────────────────────
    XA   = [ones(r,1), S*Bsm, Nuis];
    kerA = 1 + (1:size(Bsm,2));  Dk = diffmat(size(Bsm,2), 2);
    PkA  = zeros(size(XA,2));  PkA(kerA,kerA) = Dk'*Dk;
    lam_L2 = logspace(-1, 6, 12);
    [lambda_A, cvllA] = select_lambda_l2(XA, yev, PkA, lam_L2, rowfold, o.cv_folds);
    [bA, VA] = penalized_logit(XA, yev, lambda_A*PkA);
    kernel_A    = Bsm * bA(kerA);
    kernel_A_se = sqrt(sum((Bsm * VA(kerA,kerA)) .* Bsm, 2));
    [devA, dfA, aicA] = model_metrics_l2(XA, yev, bA, VA);

    % ── (B) L1 exponential ; (C) L1 flexible ──────────────────────────────────
    [~, kernel_B, kernel_B_se, sel_tau, lambda_B, cvllB, devB, dfB, aicB] = ...
        fit_l1(S, Bexp, Nuis, yev, rowfold, o.cv_folds, is_exp, tau_grid, o.sel_rule, 'B: L1 exp', o.verbose);
    [~, kernel_C, kernel_C_se, ~, lambda_C, cvllC, devC, dfC, aicC] = ...
        fit_l1(S, Bflex, Nuis, yev, rowfold, o.cv_folds, false(1,size(Bflex,2)), [], o.sel_rule, 'C: L1 flex', o.verbose);

    out = struct( ...
        't_rel',t_rel, 'caus_mask',caus_mask, 'lag_fr',lag_fr, 'fps',fps, ...
        'tau_grid',tau_grid, 'is_exp',is_exp, 'sel_tau',sel_tau, ...
        'kernel_A',kernel_A, 'kernel_A_se',kernel_A_se, ...
        'kernel_B',kernel_B, 'kernel_B_se',kernel_B_se, ...
        'kernel_C',kernel_C, 'kernel_C_se',kernel_C_se, ...
        'devN',devN, 'dfN',dfN, 'aicN',aicN, 'cv_devN',cv_devN, ...
        'devA',devA, 'dfA',dfA, 'aicA',aicA, 'cv_devA',-2*cvllA, 'lambda_A',lambda_A, ...
        'devB',devB, 'dfB',dfB, 'aicB',aicB, 'cv_devB',-2*cvllB, 'lambda_B',lambda_B, ...
        'devC',devC, 'dfC',dfC, 'aicC',aicC, 'cv_devC',-2*cvllC, 'lambda_C',lambda_C, ...
        'n_events',sum(yev), 'n_intervals',r);

    % baseline hazard b0(time-in-freeze), sm-independent (from model A): the
    % "probability of breaking at time t" separated from the sm kernel.
    base_cols = 1 + size(Bsm,2) + (1:size(Bb,2));
    out.base_logit = bA(1) + Bb * bA(base_cols);
    out.tinf       = tinf;
end

% ════════════════════════════════════════════════════════════════════════
% Local functions (mirrors of hazard_kernel_sparse.m)
% ════════════════════════════════════════════════════════════════════════
function o = default_opts(opts)
    o = struct('nb_kernel',12, 'nb_acausal',4, 'bump_width',1.25, 'nb_baseline',6, ...
               'n_exp',10, 'cv_folds',5, 'kernel_past_s',5, 'sel_rule','argmax', 'verbose',false);
    if nargin >= 1 && ~isempty(opts)
        f = fieldnames(opts);
        for i = 1:numel(f), o.(f{i}) = opts.(f{i}); end
    end
end

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

function [b, kernel, kernel_se, sel_tau, lambda, cvll, dev, df, aic, ker_cols] = ...
        fit_l1(S, Bk, Nuis, y, rowfold, K, is_exp_cols, tau_grid, sel_rule, name, verbose)
    if nargin < 9  || isempty(sel_rule), sel_rule = 'argmax'; end
    if nargin < 10, name = 'L1'; end
    if nargin < 11, verbose = false; end
    r        = size(S,1);
    X        = [ones(r,1), S*Bk, Nuis];
    nb_k     = size(Bk,2);
    ker_cols = 1 + (1:nb_k);
    penmask  = false(1, size(X,2));  penmask(ker_cols) = true;

    b_nuis   = penalized_logit(X(:, ~penmask), y, 0);
    mu0      = 1 ./ (1 + exp(-min(max(X(:,~penmask)*b_nuis,-30),30)));
    lam_max  = max(abs(X(:,penmask)' * (y - mu0)));
    lam_grid = lam_max * logspace(0, -3, 12);

    llf = nan(K, numel(lam_grid));
    for k = 1:K
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

    natom = nan(1,numel(lam_grid));  Bfull = cell(1,numel(lam_grid));  bw = [];
    for il = 1:numel(lam_grid)
        bw = penalized_logit_l1(X, y, lam_grid(il), penmask, bw);
        Bfull{il} = bw;
        natom(il) = sum(abs(bw(ker_cols)) > 1e-8);
    end

    [~, iarg] = max(cvll_lam);
    thr   = cvll_lam(iarg) - cvse_lam(iarg);
    i1se  = find(cvll_lam >= thr, 1, 'first');
    if strcmpi(sel_rule, '1se'), isel = i1se; else, isel = iarg; end
    lambda = lam_grid(isel);  cvll = cvll_lam(isel);  b = Bfull{isel};

    if verbose
        fprintf('\n── %s  λ-path  (* = argmax, s = 1-SE; sel=%s) ──\n', name, sel_rule);
        fprintf('%3s %11s %13s %7s\n', '#', 'lambda', 'CVll/obs', 'natom');
        for il = 1:numel(lam_grid)
            mk = '  '; if il==iarg, mk(1)='*'; end; if il==i1se, mk(2)='s'; end
            fprintf('%3d %11.3g %13.5f %7d %s\n', il, lam_grid(il), cvll_lam(il), natom(il), mk);
        end
    end

    bk     = b(ker_cols);
    kernel = Bk * bk;
    eta = min(max(X*b, -30), 30);  mu = min(max(1./(1+exp(-eta)),1e-9),1-1e-9);
    dev = -2 * sum(y.*log(mu) + (1-y).*log(1-mu));
    df  = sum(abs(b) > 1e-8);
    aic = dev + 2*df;

    active = ker_cols(abs(bk) > 1e-8) - 1;
    if isempty(active)
        kernel_se = zeros(size(kernel));  sel_tau = [];
    else
        Xdb = [ones(r,1), S*Bk(:,active), Nuis];
        [~, Vdb] = penalized_logit(Xdb, y, 1e-6*eye(size(Xdb,2)));
        kcol  = 1 + (1:numel(active));
        Bact  = Bk(:, active);
        kernel_se = sqrt(sum((Bact * Vdb(kcol,kcol)) .* Bact, 2));
        if ~isempty(tau_grid) && any(is_exp_cols)
            sel_tau = tau_grid(active(is_exp_cols(active)));
        else
            sel_tau = [];
        end
    end
end

function [lambda, cvll] = select_lambda_l2(X, y, Pk, lam_grid, rowfold, K)
    cvll_lam = nan(size(lam_grid));  cvse_lam = nan(size(lam_grid));
    for il = 1:numel(lam_grid)
        [cvll_lam(il), llf] = pen_grouped_cv(X, y, lam_grid(il)*Pk, rowfold, K);
        cvse_lam(il) = std(llf, 'omitnan') / sqrt(sum(~isnan(llf)));
    end
    [~, ibest] = max(cvll_lam);
    thr   = cvll_lam(ibest) - cvse_lam(ibest);
    i1se  = find(cvll_lam >= thr, 1, 'last');
    lambda = lam_grid(i1se);  cvll = cvll_lam(i1se);
end

function [dev, df, aic] = model_metrics_l2(X, y, b, V)
    eta = min(max(X*b, -30), 30);  mu = min(max(1./(1+exp(-eta)),1e-9),1-1e-9);
    dev = -2 * sum(y.*log(mu) + (1-y).*log(1-mu));
    w   = max(mu.*(1-mu), 1e-9);
    df  = trace(V * (X' * (X .* w)));
    aic = dev + 2*df;
end

function B = raised_cosine_basis(x, nb, width_mult, anchor)
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

function [b, support] = penalized_logit_l1(X, y, lambda, penmask, b0)
    [~, p] = size(X);
    if nargin < 5 || isempty(b0), b = zeros(p,1); else, b = b0(:); end
    penmask = logical(penmask(:));
    for outer = 1:50 %#ok<NASGU>
        eta = min(max(X*b, -30), 30);
        mu  = 1 ./ (1 + exp(-eta));
        w   = max(mu .* (1 - mu), 1e-6);
        z   = eta + (y - mu) ./ w;
        wxx = sum((X.^2) .* w, 1)';
        rres = z - X*b;
        b_old = b;
        for inner = 1:200 %#ok<NASGU>
            bmax = 0;
            for j = 1:p
                rho = wxx(j)*b(j) + sum(w .* X(:,j) .* rres);
                if penmask(j)
                    bj = sign(rho) * max(abs(rho) - lambda, 0) / max(wxx(j), eps);
                else
                    bj = rho / max(wxx(j), eps);
                end
                db = bj - b(j);
                if db ~= 0, rres = rres - X(:,j)*db; b(j) = bj; bmax = max(bmax, abs(db)); end
            end
            if bmax < 1e-7, break; end
        end
        if max(abs(b - b_old)) < 1e-7, break; end
    end
    support = find(abs(b) > 1e-8);
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
