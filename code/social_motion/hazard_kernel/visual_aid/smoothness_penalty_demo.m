clearvars; close all; clc;
% ════════════════════════════════════════════════════════════════════════
% smoothness_penalty_demo — VISUAL AID
%
% What the L2 smoothness penalty ||D b||^2 "sees". D is the 2nd-difference
% operator (rows [1 -2 1]); D b is the local curvature of the coefficient
% sequence, and ||D b||^2 = Σ_k (b_{k-1} - 2 b_k + b_{k+1})^2 is total
% roughness. Straight ramp -> 0 (linear trends are free); zigzag -> large.
% ════════════════════════════════════════════════════════════════════════

K = 8;
D = diffmat(K, 2);                       % (K-2) × K  second-difference operator
disp('D (2nd-difference operator):'); disp(D);

% three coefficient sequences to score
b_ramp   = linspace(0, 1, K)';                       % straight ramp  -> Db = 0
b_smooth = exp(-((1:K)' - 4.5).^2 / 5);              % smooth bump    -> small
b_zig    = mod((1:K)', 2);                           % zigzag 1 0 1 0 -> large

cases = { b_ramp,  '';
          b_smooth,'';
          b_zig,   '' };
nC = size(cases,1);

% ── print the penalty for each ───────────────────────────────────────────────
fprintf('\n%-14s  ||Db||^2\n', 'sequence');
for c = 1:nC
    b = cases{c,1};  Db = D*b;
    fprintf('%-14s  %8.3f     Db = [%s]\n', cases{c,2}, Db'*Db, ...
        strtrim(sprintf('%5.2f ', Db)));
end

% ── figure: top row = coefficients b, bottom row = second differences Db ─────
kb  = (1:K)';           % coefficient index
kd  = (2:K-1)';         % interior index where Db lives
cB  = [0.25 0.45 0.75];
cD  = [0.85 0.37 0.12];

fh = figure('Color','w','Position',[60 60 900 400]);
tl = tiledlayout(2, nC, 'TileSpacing','compact', 'Padding','compact'); %#ok<NASGU>

for c = 1:nC
    b = cases{c,1};  Db = D*b;  pen = Db'*Db;

    % top: the coefficient sequence b_k
    ax = nexttile(c); hold(ax,'on');
    stem(ax, kb, b, 'filled', 'Color', cB, 'MarkerFaceColor', cB, 'LineWidth', 1.5);
    plot(ax, kb, b, ':', 'Color', cB);
    ylim(ax, [min(-0.2, min(b)-0.2), max(b)+0.3]); xlim(ax,[0.5 K+0.5]);
    title(ax, sprintf('%s||Db||^2 = %.2f ', cases{c,2}, pen), 'FontWeight','normal');
    if c==1, ylabel(ax, 'coefficient  b_k'); end
    set(ax,'Box','off','FontSize',24,'XTick',1:K,'TickDir','out','Layer','top');
    apply_generic(gca)

    % bottom: the second differences Db_k = b_{k-1} - 2 b_k + b_{k+1}
    ax2 = nexttile(nC + c); hold(ax2,'on');
    stem(ax2, kd, Db, 'filled', 'Color', cD, 'MarkerFaceColor', cD, 'LineWidth', 1.5);
    plot(ax2, [0.5 K+0.5], [0 0], 'k-', 'LineWidth', 0.6);
    yl = max(1, max(abs(Db))*1.3); ylim(ax2, [-yl yl]); xlim(ax2,[0.5 K+0.5]);
    if c==1, ylabel(ax2, 'curvature  (D b)_k'); end
    xlabel(ax2, 'coefficient index k');
    set(ax2,'Box','off','FontSize',24,'XTick',1:K,'TickDir','out','Layer','top');
    apply_generic(gca)
end

% ── save (PNG preview + vector PDF via the repo helper) ──────────────────────
outdir = fileparts(mfilename('fullpath'));
exportgraphics(fh, fullfile(outdir,'smoothness_penalty_demo.png'), 'Resolution',200);
paths = path_generator('folder','social_motion/hazard_kernel/visual_aid','imfirst',false);
exporter(fh, paths, 'smoothness_penalty_demo.pdf');   % vector PDF -> paths.fig
fprintf('\nsaved: %s\n       %s\n', fullfile(outdir,'smoothness_penalty_demo.png'), ...
    fullfile(paths.fig,'smoothness_penalty_demo.pdf'));

% ── local: 2nd-difference operator (same as the fit scripts) ─────────────────
function D = diffmat(n, ord)
    D = eye(n);
    for k = 1:ord, D = diff(D); end
end
