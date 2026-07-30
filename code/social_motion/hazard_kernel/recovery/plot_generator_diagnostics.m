function [fh_dur, fh_sm] = plot_generator_diagnostics(bl_all, bl_t, motion_cache, beta_list, cfg, P)
%PLOT_GENERATOR_DIAGNOSTICS  Two sanity figures for the DDM recovery generator.
%
%   [fh_dur, fh_sm] = plot_generator_diagnostics(bl_all, bl_t, motion_cache, ...
%                                                beta_list, cfg, P)
%
%   bl_all : {nBeta x 1} cell of SYNTHETIC bout tables from
%            generate_accumulator_dataset (needs .fly .onsets .dur_frames
%            .is_censored), one per entry of beta_list.
%   bl_t   : the REAL template table (uses .durations as the reference).
%   cfg    : the build_hazard_design config (uses .fps .lag_fr .sm_mu .sm_sd).
%   P      : the generator struct (uses .mu0 .theta .n_aug for annotation).
%
% FIGURE 1 (fh_dur) — simulated freeze-duration distribution per beta, against the
%   real distribution. These vary with BETA only: mask_preonset / nb_baseline are
%   design-side options and never touch the simulation, so re-running the fit with
%   different estimator settings reproduces these panels exactly.
%
% FIGURE 2 (fh_sm) — grand-average social motion aligned to freeze OFFSET, drawn
%   three ways that correspond to the three ways a design matrix can treat lags
%   reaching back past freeze onset:
%     'in-freeze only'  NaN outside the freeze, nanmean  -> the honest average,
%                       but computed over a SHRINKING subset as you go back:
%                       only long bouts have 3 s of in-freeze history.
%     'pre-onset -> 0'  zero-filled outside (mask_preonset = true). Because the
%                       fill is applied on exactly the rows that are too young,
%                       the column doubles as a bout-age indicator.
%     'all history'     no masking (mask_preonset = false): real sm from before
%                       onset, and from outside the loom period entirely.
%   Everything is in the estimator's z-scored units, so the zero line IS the
%   mask fill value. The dashed reference marks raw sm = 0 -- note it is NOT at
%   z = 0, i.e. "masked" means "average social motion", not "no social motion".

    fps      = cfg.fps;
    nBeta    = numel(beta_list);
    back_fr  = max(cfg.lag_fr);              % causal window of the design
    post_fr  = round(0.5 * fps);             % a little context after the un-freeze
    col_bg   = [0.55 0.55 0.55];             % all history
    col_zero = [0.85 0.40 0.20];             % pre-onset -> 0
    col_frz  = [0.20 0.50 0.80];             % in-freeze only

    % ════════════════════════════════════════════════════════════════════════
    % FIGURE 1 — duration distributions
    % ════════════════════════════════════════════════════════════════════════
    real_s = double(bl_t.durations) / fps;
    all_s  = real_s;
    for m = 1:nBeta, all_s = [all_s; bl_all{m}.dur_frames / fps]; end %#ok<AGROW>
    xmax   = prctile(all_s, 99);
    edges  = linspace(0, xmax, 60);

    fh_dur = figure('Color','w','Position',[60 60 330*nBeta 360]);
    tiledlayout(1, nBeta, 'TileSpacing','compact','Padding','compact');
    for m = 1:nBeta
        bl  = bl_all{m};
        d   = bl.dur_frames / fps;
        unc = ~bl.is_censored;
        nexttile; hold on

        histogram(d(unc), edges, 'Normalization','pdf', ...
            'FaceColor',[0.35 0.55 0.75], 'EdgeColor','none', 'FaceAlpha', 0.85);
        xline(median(d(unc)), '-',  'Color',[0.35 0.55 0.75], 'LineWidth',1.6);

        xlim([0 xmax]);
        xlabel('Freeze Duration (s)');
        if m == 1
            ylabel('pdf');
        end
        apply_generic(gca, 'xlim', [0 10.5], 'ylim', [0 1])
    end

    % ════════════════════════════════════════════════════════════════════════
    % FIGURE 2 — social motion aligned to freeze offset, three masking rules
    % ════════════════════════════════════════════════════════════════════════
    w    = (-back_fr : post_fr)';            % frames relative to un-freeze
    t_w  = w / fps;
    z0   = -cfg.sm_mu / cfg.sm_sd;           % where raw sm = 0 sits in z units

    % ── pass 1: accumulate the three averages per beta ──────────────────────
    mu_all = nan(numel(w), nBeta);           % no masking
    mu_zer = nan(numel(w), nBeta);           % pre-onset -> 0
    mu_frz = nan(numel(w), nBeta);           % in-freeze only (nanmean)
    ncon   = nan(numel(w), nBeta);           % bouts contributing to mu_frz
    nbout  = nan(1, nBeta);
    for m = 1:nBeta
        bl   = bl_all{m};
        keep = ~bl.is_censored;              % offset only defined for real un-freezes
        fly  = bl.fly(keep);  on = bl.onsets(keep);
        off  = on + bl.dur_frames(keep) - 1;
        n    = numel(off);  nbout(m) = n;

        A = nan(n, numel(w));                % z-scored sm, NaN outside the recording
        M = false(n, numel(w));              % true where the sample is INSIDE the freeze
        for i = 1:n
            sm  = motion_cache(fly(i));  sm = sm(:);  L = numel(sm);
            idx = off(i) + w;
            ok  = idx >= 1 & idx <= L;
            A(i,ok) = (sm(idx(ok)) - cfg.sm_mu) / cfg.sm_sd;
            M(i,:)  = ok & idx >= on(i) & idx <= off(i);
        end

        Afrz = A;  Afrz(~M) = NaN;
        Azer = A;  Azer(~M) = 0;
        mu_all(:,m) = mean(A,    1, 'omitnan')';
        mu_zer(:,m) = mean(Azer, 1, 'omitnan')';
        mu_frz(:,m) = mean(Afrz, 1, 'omitnan')';
        ncon(:,m)   = sum(~isnan(Afrz), 1)';
    end

    % common scales so the panels are directly comparable
    lo = min([mu_all(:); mu_zer(:); mu_frz(:); z0], [], 'omitnan');
    hi = max([mu_all(:); mu_zer(:); mu_frz(:); 0 ], [], 'omitnan');
    pad = 0.08 * (hi - lo);
    ylimL = [lo - pad, hi + pad];
    ylimR = [0, max(nbout) * 1.05];

    % ── pass 2: draw ────────────────────────────────────────────────────────
    fh_sm = figure('Color','w','Position',[60 60 450*nBeta 420]);
    tl = tiledlayout(1, nBeta, 'TileSpacing','compact','Padding','compact');
    for m = 1:nBeta
        nexttile; hold on
        apply_generic(gca, 'xlim', [-3 0])
        fill([0 -0.5 -0.5 0], [2 2 -2 -2], 'k', 'FaceAlpha', 0.1)
        yyaxis right                          % draw first so it sits behind
        plot(t_w, ncon(:,m), '--', 'Color',[0.2 0.2 0.2], 'LineWidth',1.2);
        ylim(ylimR);  set(gca,'YColor',[0.2 0.2 0.2]);
        % yyaxis axes have a 2x1 YAxis (left+right ruler), so index the side --
        % axR.YAxis.LineWidth = ... errors on the non-scalar array.
        axR = gca; axR.YAxis(2).LineWidth = 2.5;
        axR.FontSize = 24;

        if m == nBeta, ylabel('# Bouts'); end

        yyaxis left
        g(1) = plot(t_w, mu_all(:,m), '-', 'Color',col_bg,   'LineWidth', 2.2);
        g(2) = plot(t_w, mu_zer(:,m), '-', 'Color',col_zero, 'LineWidth', 2.2);
        g(3) = plot(t_w, mu_frz(:,m), '-', 'Color',col_frz,  'LineWidth', 2.2);
        yline(0,  'k-', 'LineWidth', 2);
        xline(0,  'k-', 'LineWidth', 2);

        ylim(ylimL);  set(gca,'YColor','k'); axL = gca; axL.YAxis(1).LineWidth = 2.5;
        axL.FontSize = 24;
        if m == 1, ylabel('Social Motion (z)'); end

        xlim([min(t_w) max(t_w)]);
        xlabel('Time to unfreeze (s)');
        title(sprintf('\\beta_{gen}=%.3g', beta_list(m)), 'FontSize', 18);

    end
    lg = legend(g, {'all history (mask off)', 'pre-onset = 0 (mask on)', ...
                    'in-freeze only (nanmean)'}, 'Orientation','horizontal','box','off');
    lg.Layout.Tile = 'south'; lg.FontSize = 20;

end
