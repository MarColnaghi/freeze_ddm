function f = memory_feature(D, motion_cache, spec)
% MEMORY_FEATURE  One causal memory feature of social motion, aligned to the
% at-risk intervals in D.  Returns a z-scored column f (D.n_int x 1).
%
%   spec.type:
%     'inst'  — instantaneous m(t)
%     'box'   — trailing boxcar mean over spec.timescale seconds  M_T(t)
%     'exp'   — leaky-integrator state (causal EWMA), leak τ = spec.timescale s:
%                 x[t] = a·x[t-1] + (1-a)·m[t],  a = exp(-1/(τ·fps))
%               This is exactly the leaky-accumulator memory; best τ = the leak.
%     'slope' — causal linear-regression slope of m over spec.timescale seconds
%               (derivative / "is motion increasing?"). Sign: + = rising toward now.
%
% Each fly's FULL trace is z-scored (global sm_mu/sm_sd from D), missing→0
% (= baseline), then filtered causally — so the memory has equilibrated by the
% time a freeze starts, and integrates pre-freeze history too (as a real leaky
% integrator would). The filtered trace is sampled at each interval's frame.
%
% NOTE: CV-LL of a single-predictor logistic model is invariant to linear
% rescaling of the predictor, so the z-scoring is only for interpretability —
% boxcar-vs-exp comparisons are scale-free and depend only on the feature SHAPE.

    fps   = D.fps;
    flies = unique(D.fly);
    cache = containers.Map('KeyType', 'double', 'ValueType', 'any');

    for i = 1:numel(flies)
        fl = flies(i);
        m  = motion_cache(fl); m = m(:);
        m  = (m - D.sm_mu) / D.sm_sd;     % global z-score (preserves cross-fly scale)
        m(isnan(m)) = 0;                  % missing sm → baseline (no social drive)

        switch spec.type
            case 'inst'
                g = m;
            case 'box'
                Tfr = max(1, round(spec.timescale * fps));
                g = movmean(m, [Tfr - 1, 0]);          % trailing (causal) mean
            case 'exp'
                a = exp(-1 / (spec.timescale * fps));   % EWMA = leaky integrator
                g = filter(1 - a, [1, -a], m);
            case 'slope'
                g = causal_slope(m, max(2, round(spec.timescale * fps)));
            otherwise
                error('memory_feature:type', 'unknown spec.type ''%s''', spec.type);
        end
        cache(fl) = g;
    end

    f = nan(D.n_int, 1);
    for r = 1:D.n_int
        g = cache(D.fly(r)); t = D.tabs(r);
        if t >= 1 && t <= numel(g), f(r) = g(t); end
    end
    f = (f - mean(f, 'omitnan')) / std(f, 'omitnan');   % z-score (interpretability)
end

function s = causal_slope(m, Tfr)
% Slope of an ordinary least-squares line fit to the last Tfr samples, evaluated
% causally at each t (FIR filter with the OLS slope weights). Edge (<Tfr) → NaN.
    k  = (0:Tfr-1)';            % lag into the past
    tt = -k;                    % time axis (now = 0, past < 0)
    w  = (tt - mean(tt));
    w  = w / sum(w.^2);         % OLS slope weights: slope = Σ w_k · m(t-k)
    s  = filter(w, 1, m);       % w(1) multiplies m(t), w(2) m(t-1), ...
    s(1:Tfr-1) = NaN;
end
