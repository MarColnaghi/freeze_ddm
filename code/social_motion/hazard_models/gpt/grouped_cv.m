function [ll, ll_folds] = grouped_cv(T, formula, rowfold, K)
% GROUPED_CV  Mean held-out Bernoulli log-likelihood/obs, folds grouped by bout.
%   Returns the across-fold mean `ll` and the per-fold means `ll_folds` (Kx1).
%   Per-fold means (not pooled) so a paired SE across folds is well defined.
    ll_folds = nan(K, 1);
    for k = 1:K
        te = rowfold == k; tr = ~te;
        if ~any(te) || ~any(tr), continue; end
        mm = fitglm(T(tr, :), formula, 'Distribution', 'binomial');
        p  = predict(mm, T(te, :));
        p  = min(max(p, 1e-9), 1 - 1e-9);
        y  = T.y(te);
        ll_folds(k) = mean(y .* log(p) + (1 - y) .* log(1 - p));
    end
    ll = mean(ll_folds, 'omitnan');
end
