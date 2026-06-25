function out = fit_cv(D, F, names)
% FIT_CV  Discrete-time logistic hazard fit with a shared baseline + controls and
% an optional set of memory features F (D.n_int x p, names = cellstr).
%
%   logit P(unfreeze) = baseline(time-in-freeze) + controls + Σ F
%
% Returns: out.aic, out.cvll (grouped CV-LL/obs), out.cvll_folds (per fold),
%          out.insample_ll (full-data LL/obs), out.k (n params), out.mdl.

    if nargin < 2, F = []; end
    if nargin < 3, names = {}; end

    T = D.ctrl_tbl; T.y = D.yev;
    nb = size(D.Bb, 2);
    for j = 1:nb, T.(sprintf('h%d', j)) = D.Bb(:, j); end
    ht = strjoin(arrayfun(@(j) sprintf('h%d', j), 1:nb, 'uni', 0), ' + ');

    fterms = '';
    for j = 1:size(F, 2)
        T.(names{j}) = F(:, j);
        fterms = [fterms ' + ' names{j}]; %#ok<AGROW>
    end

    rhs = ['1 + ' ht D.ctrl_terms fterms];
    mdl = fitglm(T, ['y ~ ' rhs], 'Distribution', 'binomial');

    out.mdl = mdl;
    out.aic = mdl.ModelCriterion.AIC;
    out.k   = mdl.NumEstimatedCoefficients;
    [out.cvll, out.cvll_folds] = grouped_cv(T, ['y ~ ' rhs], D.rowfold, D.cv_folds);

    p = predict(mdl, T); p = min(max(p, 1e-9), 1 - 1e-9);
    out.insample_ll = mean(T.y .* log(p) + (1 - T.y) .* log(1 - p));
end
