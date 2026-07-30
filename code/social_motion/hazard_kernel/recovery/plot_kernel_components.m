function fh = plot_kernel_components(src, opts, lag_fr, fps, beta_list)
%PLOT_KERNEL_COMPONENTS  The weighted basis components b_j*B(:,j) behind kernels A and B.
%
%   fh = plot_kernel_components(R, opts, lag_fr, fps, beta_list)
%   fh = plot_kernel_components('results/ddm_kernel_recovery_20260730_095354.mat')
%
%   Same idea as the third panel of visual_aid/kernel_fit_single.m (model A,
%   raised cosines) and visual_aid/kernel_fit_exp.m (model B, exponential bank):
%   thin coloured lines are the individual weighted basis functions, the thick
%   black line is their sum -- the kernel ddm_kernel_recovery.m plots.
%
%   Rows: A = raised cosines + L2 curvature ; B = exponential bank + L1.
%   Cols: one per entry of beta_list.
%
%   fit_kernels_sparse returns the kernels but not the coefficients, so b is
%   recovered by least squares from kernel = B*b. That is EXACT rather than an
%   approximation -- the kernel lies in the column span of B by construction. For
%   B the solve is restricted to the atoms the lasso actually kept (res.sel_tau):
%   a min-norm solve on a bank that is 99.8% collinear would smear one sparse
%   kernel across every atom and misrepresent what L1 did.

    % ── accept either the in-memory sweep or a saved run ────────────────────
    if ischar(src) || isstring(src)
        S = load(char(src), 'R', 'opts', 'lag_fr', 'beta_list');
        R = S.R;  opts = S.opts;  lag_fr = S.lag_fr;  beta_list = S.beta_list;
        fps = R(1).res.fps;
    else
        R = src;
    end
    nBeta = numel(R);

    % ── rebuild the two bases (verbatim from fit_kernels_sparse) ────────────
    Bsm  = raised_cosine_basis(lag_fr, opts.nb_kernel + opts.nb_acausal, opts.bump_width);
    % runs saved before tau_lo existed used the old hard-coded 0.05 s lower bound
    if isfield(opts,'tau_lo'), tau_lo = opts.tau_lo; else, tau_lo = 0.05; end
    tau_grid_s = logspace(log10(tau_lo), log10(opts.kernel_past_s), opts.n_exp);
    [Bexp, tau_grid, is_exp] = exp_basis(lag_fr, tau_grid_s, fps, opts.nb_acausal, opts.bump_width);

    KA = size(Bsm,2);  KB = size(Bexp,2);
    cmapA  = cbrewer2('Spectral', KA);   % ordered by bump centre / by tau
    cmapB  = cbrewer2('Spectral', KB);

    [t_ax, ord] = sort(-lag_fr / fps);

    % ── pass 1: recover coefficients, collect components, find common scales ─
    compA = cell(nBeta,1);  compB = cell(nBeta,1);
    bA    = cell(nBeta,1);  bB    = cell(nBeta,1);
    kerA  = zeros(numel(lag_fr), nBeta);  kerB = kerA;
    for m = 1:nBeta
        rr = R(m).res;
        kerA(:,m) = rr.kernel_A;  kerB(:,m) = rr.kernel_B;

        bA{m}    = Bsm \ rr.kernel_A;              % full rank: exact
        compA{m} = Bsm .* bA{m}(:).';

        % restrict the B solve to the L1 support
        act = false(1, KB);
        if ~isempty(rr.sel_tau)
            act(is_exp) = ismembertol(tau_grid, rr.sel_tau(:).', 1e-9);
        end
        if opts.nb_acausal > 0, act(~is_exp) = true; end   % acausal cols are unpenalised
        bb = zeros(KB,1);
        if any(act), bb(act) = Bexp(:,act) \ rr.kernel_B; end
        bB{m}    = bb;
        compB{m} = Bexp .* bb(:).';
    end
    ylA = pad_lim([compA{:}, kerA]);
    ylB = pad_lim([compB{:}, kerB]);

    % ── pass 2: draw ───────────────────────────────────────────────────────
    fh = figure('Color','w','Position',[60 60 430*nBeta 500]);
    tiledlayout(2, nBeta, 'TileSpacing','compact','Padding','compact');

    % Every call is given its axes handle explicitly: with a 2-row tiledlayout,
    % relying on gca puts the annotations in whichever tile was touched last.
    for m = 1:nBeta                                   % ── row 1: model A ──
        ax = nexttile(m); hold(ax,'on');
        yline(ax, 0, 'k-', 'LineWidth',1.5);
        for j = 1:KA
            plot(ax, t_ax, compA{m}(ord,j), '-', 'Color', cmapA(j,:), 'LineWidth', 1.2);
        end
        plot(ax, t_ax, kerA(ord,m), 'k-', 'LineWidth', 2.6);
        xlim(ax, [min(t_ax) max(t_ax)]);  ylim(ax, ylA);
        title(ax, sprintf('\\beta_{gen}=%.3g', beta_list(m)));
        if m == 1, ylabel(ax, 'A-L2'); end
        apply_generic(ax, 'ylim', [0 0.2]);
        inset_bars(fh, ax, bA{m}, cmapA, 'bump  (past \rightarrow 0)');
    end

    for m = 1:nBeta                                   % ── row 2: model B ──
        ax = nexttile(nBeta + m); hold(ax,'on');
        yline(ax, 0, 'k-', 'LineWidth',1.5);
        for j = find(abs(bB{m}) > 0).'
            plot(ax, t_ax, compB{m}(ord,j), '-', 'Color', cmapB(j,:), 'LineWidth', 1.6);
        end
        plot(ax, t_ax, kerB(ord,m), 'k-', 'LineWidth', 2.6);
        xlim(ax, [min(t_ax) max(t_ax)]);  ylim(ax, ylB);
        % which atoms survived L1 -- the whole point of the sparse fit
        if isempty(R(m).res.sel_tau)
            txt = 'no atoms selected';
        else
            txt = sprintf('\\tau = %s s', num2str(sort(R(m).res.sel_tau), '%.2f  '));
        end
        % DATA units, not normalized: exportgraphics re-lays out the tiledlayout
        % and normalized-unit text does not follow its axes through that.
        text(ax, max(t_ax) - 0.02*range(t_ax), ylB(2) - 0.06*diff(ylB), txt, ...
            'HorizontalAlignment','right', 'VerticalAlignment','top', 'FontSize', 14);
        xlabel(ax, 'Time to unfreeze (s)');
        if m == 1, ylabel(ax, 'B-L1'); end
        apply_generic(ax, 'ylim', [0 0.2]);
        inset_bars(fh, ax, bB{m}, cmapB, '\tau  (fast \rightarrow slow)');
    end
end

% ════════════════════════════════════════════════════════════════════════
function inset_bars(fh, ax, b, cmap, xlab)
% Small bar panel of the basis coefficients b_j, in the top-left of AX (where the
% kernel is flat -- it rises towards t = 0 on the right). Bars are coloured to
% match the component curves, so a bar and its curve are the same object.
%
% Placed from AX.Position AFTER apply_generic has run, since that can change the
% axes extent. Every bar is drawn, including the zeros: for the L1 fits the gaps
% ARE the result, and dropping them would hide how sparse the fit is.
    K   = numel(b);
    p   = ax.Position;
    pos = [p(1) + 0.10*p(3), p(2) + 0.56*p(4), 0.40*p(3), 0.36*p(4)];
    axi = axes('Parent', fh, 'Position', pos);  hold(axi,'on');

    for j = 1:K
        bar(axi, j, b(j), 0.8, 'FaceColor', cmap(j,:), 'EdgeColor','none');
    end
    plot(axi, [0.4 K+0.6], [0 0], 'k-', 'LineWidth', 0.8);

    xlim(axi, [0.4 K+0.6]);
    hi = max(abs(b));  if hi <= 0, hi = 1; end
    ylim(axi, [min(0, 1.15*min(b)), 1.15*hi]);
    % label the peak rather than rounding to 2 dp -- at beta = 0 the coefficients
    % are ~1e-3, which rounds to 0 and gives a non-increasing tick vector
    set(axi, 'XTick', [], 'YTick', [0 hi], 'YTickLabel', {'0', sprintf('%.2g', hi)}, ...
        'Box','off', 'Color','none', 'FontSize', 11, 'TickDir','out', 'LineWidth', 1);
    xlabel(axi, xlab, 'FontSize', 11);
    ylabel(axi, 'b_j', 'FontSize', 11);
end

% ════════════════════════════════════════════════════════════════════════
function yl = pad_lim(V)
    lo = min(V(:));  hi = max(V(:));  pad = 0.06 * max(hi - lo, eps);
    yl = [lo - pad, hi + pad];
end

% ── local basis builders: verbatim copies from fit_kernels_sparse.m, which
% keeps them private. They must stay in step with that file.
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
