function D = build_hazard_design(bl, motion_cache, opts)
% BUILD_HAZARD_DESIGN  Person-period (discrete-time) hazard design for freeze
% termination, shared by every step of the memory-timescale pipeline.
%
%   D = build_hazard_design(bl, motion_cache, opts)
%
% Each at-risk interval (sampled every dt_frames from the left-truncation entry
% time to the bout end) becomes one row. We DO NOT build a lagged kernel design
% here — instead we store the absolute frame index of each interval (D.tabs) so
% memory_feature.m can sample any causal filter of the fly's sm trace at that
% moment. That keeps every memory model (instantaneous / boxcar / exponential) on
% the identical rows, events, controls, and CV folds — so all ΔCV-LL are paired.
%
% bl (table), one row per freeze bout, must contain:
%   fly, onsets, dur_frames, is_censored, moving_flies, sloom, bout_id
% motion_cache : containers.Map, fly id -> per-frame sum_motion vector.
%
% opts (optional): fps(60) dt_frames(6) nb_baseline(6) cv_folds(5)
%                  trunc_point(30) rng_seed(7)
%
% D fields: tabs, tinf(s), yev, gid, fly, mov, slo, Bb (baseline basis, intercept
%   column dropped), ctrl_tbl, ctrl_terms, rowfold, sm_mu, sm_sd, fps, dt_frames,
%   cv_folds, n_int, n_events.

    fps         = geto(opts, 'fps', 60);
    dt_frames   = geto(opts, 'dt_frames', 6);
    nb_baseline = geto(opts, 'nb_baseline', 6);
    cv_folds    = geto(opts, 'cv_folds', 5);
    trunc_point = geto(opts, 'trunc_point', 30);
    grid_anchor = geto(opts, 'grid_anchor', 'entry');  % 'entry' from onset | 'offset' from break
    rng_seed    = geto(opts, 'rng_seed', 7);
    rng(rng_seed);
    entry_fr = trunc_point;

    % ── global z-scoring scale for sm (one scale across flies) ───────────────
    flies = unique(bl.fly);
    nn = 0; s1 = 0; s2 = 0;
    for i = 1:numel(flies)
        v = motion_cache(flies(i)); v = v(:); v = v(~isnan(v));
        nn = nn + numel(v); s1 = s1 + sum(v); s2 = s2 + sum(v.^2);
    end
    sm_mu = s1 / nn;  sm_sd = sqrt(s2 / nn - sm_mu^2);

    % ── person-period expansion with left truncation / delayed entry ─────────
    max_rows = sum(ceil(bl.dur_frames / dt_frames)) + height(bl);
    tabs = nan(max_rows, 1); tinf = nan(max_rows, 1); yev = zeros(max_rows, 1);
    gid  = nan(max_rows, 1); flyv = nan(max_rows, 1);
    mov  = nan(max_rows, 1); slo  = nan(max_rows, 1);
    r = 0;
    for b = 1:height(bl)
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
            f = grid(gi); r = r + 1;
            tabs(r) = on + f - 1;       % absolute frame of this interval
            tinf(r) = f / fps;          % time-in-freeze (s) for the baseline hazard
            gid(r)  = bl.bout_id(b);
            flyv(r) = bl.fly(b);
            mov(r)  = bl.moving_flies(b);
            slo(r)  = bl.sloom(b);
            % event only on the final interval of an UNCENSORED bout
            if gi == numel(grid) && ~bl.is_censored(b), yev(r) = 1; end
        end
    end
    ix = 1:r;
    D.tabs = tabs(ix); D.tinf = tinf(ix); D.yev = yev(ix);
    D.gid  = gid(ix);  D.fly  = flyv(ix); D.mov = mov(ix); D.slo = slo(ix);

    % ── baseline-hazard basis on time-in-freeze (rank-scale knots) ───────────
    % Equal-mass knots (empirical CDF of tinf); drop one bump (raised cosines sum
    % to ≈1, would duplicate the intercept fitglm adds).
    ut = tiedrank(D.tinf) / numel(D.tinf);
    Bb = raised_cosine_basis(ut, nb_baseline);
    D.Bb = Bb(:, 2:end);

    % ── controls: include mov/slo only if they vary ─────────────────────────
    ctrl_tbl = table(); ctrl_terms = '';
    if numel(unique(D.mov)) > 1, ctrl_tbl.mov = categorical(D.mov); ctrl_terms = [ctrl_terms ' + mov']; end
    if numel(unique(D.slo)) > 1, ctrl_tbl.slo = categorical(D.slo); ctrl_terms = [ctrl_terms ' + slo']; end
    D.ctrl_tbl = ctrl_tbl; D.ctrl_terms = ctrl_terms;

    % ── grouped (by bout) CV folds — shared across the whole pipeline ────────
    fold = mod(randperm(max(D.gid)), cv_folds) + 1;
    D.rowfold = fold(D.gid);

    D.fps = fps; D.dt_frames = dt_frames; D.cv_folds = cv_folds;
    D.sm_mu = sm_mu; D.sm_sd = sm_sd;
    D.n_int = r; D.n_events = sum(D.yev);
    fprintf('Design: %d at-risk intervals, %d un-freeze events (hazard %.3f), %d bouts\n', ...
        D.n_int, D.n_events, mean(D.yev), max(D.gid));
end
