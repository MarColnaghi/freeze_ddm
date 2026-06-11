function P = profile_timescale(D, motion_cache, type, grid)
% PROFILE_TIMESCALE  Fit a single-predictor memory model at each timescale in
% `grid` (seconds) and return the predictive-performance profile.
%
%   type : 'box' or 'exp' (or 'inst'/'slope' with a 1-element grid).
%
% P fields (per grid point): cvll, cvse (SE across folds), aic, insample_ll,
%   plus P.grid, P.type, and the argmax-by-CV summary P.best_idx/best_ts/best_cvll.
%
% Each model is ONE z-scored predictor on top of the shared baseline+controls, so
% the profile is a clean, multicollinearity-free, scale-free comparison of memory
% SHAPES. The location of the peak = the preferred averaging/leak timescale; the
% flatness of the peak = how well that timescale is identified.

    n = numel(grid);
    P.type = type; P.grid = grid(:);
    P.cvll = nan(n,1); P.cvse = nan(n,1); P.aic = nan(n,1); P.insample_ll = nan(n,1);
    for i = 1:n
        f   = memory_feature(D, motion_cache, struct('type', type, 'timescale', grid(i)));
        out = fit_cv(D, f, {'mem'});
        P.cvll(i)        = out.cvll;
        P.cvse(i)        = std(out.cvll_folds, 'omitnan') / sqrt(sum(~isnan(out.cvll_folds)));
        P.aic(i)         = out.aic;
        P.insample_ll(i) = out.insample_ll;
    end
    [P.best_cvll, P.best_idx] = max(P.cvll);
    P.best_ts = grid(P.best_idx);
end
