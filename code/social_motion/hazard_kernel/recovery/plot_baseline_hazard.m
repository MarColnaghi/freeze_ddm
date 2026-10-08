function fh = plot_baseline_hazard(src, bl_all, beta_list, cfg)
%PLOT_BASELINE_HAZARD  The fitted duration dependence b0(time-in-freeze), against
%the accumulator's true first-passage hazard.
%
%   fh = plot_baseline_hazard(R, bl_all, beta_list, cfg)
%   fh = plot_baseline_hazard('results/ddm_kernel_recovery_20260730_132626.mat')
%
%   One panel per beta. Solid line + band = b0(t) from model A with its 95% CI,
%   converted from per-interval log-odds to a rate per second so it is readable
%   independently of dt_frames. Points = the DDM's own hazard, measured directly
%   from the simulated durations in bl_all (events per second at risk, in
%   equal-mass duration bins) -- i.e. the truth b0(t) is trying to recover.
%
%   Returns [] if the run had use_nuisance = false: there is then no baseline
%   block in the design and base_logit is just the intercept repeated, which
%   would draw a flat line that looks like a finding rather than an absence.
%
%   READ THE TWO CURVES AS SHAPES, NOT LEVELS. b0(t) is the hazard when the kernel
%   term is zero (z-scored sm history at its mean), while the empirical hazard
%   marginalises over the sm actually seen. With a non-zero kernel those differ by
%   a vertical offset, so only the beta = 0 panel is a like-for-like calibration
%   check -- there the two should coincide.

    % ── accept either the in-memory sweep or a saved run ────────────────────
    if ischar(src) || isstring(src)
        S = load(char(src), 'R', 'bl_all', 'beta_list', 'cfg');
        R = S.R;  bl_all = S.bl_all;  beta_list = S.beta_list;  cfg = S.cfg;
    else
        R = src;
    end

    if ~isfield(R(1).res, 'base_se')
        warning(['plot_baseline_hazard: this run has no baseline block ' ...
                 '(use_nuisance = false) -- nothing to draw.']);
        fh = [];  return
    end

    nBeta = numel(R);
    dt    = cfg.dt_frames / cfg.fps;      % seconds per at-risk interval
    fps   = cfg.fps;
    col_f = [0.20 0.50 0.80];             % fitted b0
    col_e = [0.25 0.25 0.25];             % simulated truth

    fh = figure('Color','w','Position',[60 60 430*nBeta 380]);
    tiledlayout(1, nBeta, 'TileSpacing','compact','Padding','compact');

    for m = 1:nBeta
        rr = R(m).res;

        % ── fitted b0(t): one value per DISTINCT time-in-freeze, not per row ──
        % b0 is a function of tinf alone, so the ~500k design rows collapse to a
        % few hundred grid times.
        [t, ia] = unique(rr.tinf);
        lo  = rr.base_logit(ia);   se = rr.base_se(ia);
        rate    = logit2rate(lo,             dt);
        rate_hi = logit2rate(lo + 1.96*se,   dt);
        rate_lo = logit2rate(lo - 1.96*se,   dt);

        % ── the truth: discrete hazard of the simulated durations ────────────
        % Equal-mass bins. at_risk counts every bout still frozen at the bin's
        % start (censored ones included -- they are at risk, they just never
        % produce an event), so this is the standard life-table estimator.
        bl  = bl_all{m};
        d   = double(bl.dur_frames) / fps;
        unc = ~bl.is_censored;
        q   = unique(prctile(d, 0:5:100));
        ar  = arrayfun(@(k) sum(d >= q(k)),                          1:numel(q)-1);
        ev  = arrayfun(@(k) sum(unc & d >= q(k) & d < q(k+1)),       1:numel(q)-1);
        rate_emp = ev ./ max(ar,1) ./ diff(q);
        tmid     = 0.5 * (q(1:end-1) + q(2:end));

        ax = nexttile(m); hold(ax,'on');
        patch(ax, [t; flipud(t)], [rate_hi; flipud(rate_lo)], col_f, ...
            'FaceAlpha', 0.20, 'EdgeColor','none');
        h1 = plot(ax, t, rate, '-', 'Color', col_f, 'LineWidth', 2.6);
        h2 = plot(ax, tmid, rate_emp, 'o', 'Color', col_e, 'MarkerSize', 5, ...
            'MarkerFaceColor','w', 'LineWidth', 1.2);

        % per-panel y scale on purpose: the sweep is NOT rate-matched (beta lifts
        % the mean drift too), so a shared axis would flatten the beta = 0 panel
        xlim(ax, [0, prctile(d, 99)]);
        xlabel(ax, 'Time in freeze (s)');
        if m == 1
            ylabel(ax, 'Hazard (s^{-1})');
            legend(ax, [h1 h2], {'fitted b_0(t) \pm 95% CI', 'simulated DDM'}, ...
                'Location','northwest', 'box','off');
        end
        title(ax, sprintf('\\beta_{gen}=%.3g', beta_list(m)));
        apply_generic(ax,'ylim', [0 4], 'xlim', [0 5]);
    end
end

% ════════════════════════════════════════════════════════════════════════
function rate = logit2rate(lg, dt)
% Per-interval log-odds -> continuous-time rate per second. The per-interval
% probability alone would be silently tied to dt_frames.
    p    = 1 ./ (1 + exp(-min(max(lg, -30), 30)));
    rate = -log(max(1 - p, eps)) / dt;
end
