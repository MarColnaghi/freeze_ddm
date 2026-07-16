clearvars; close all; clc;
% ════════════════════════════════════════════════════════════════════════
% Visualise the two kernel bases used in hazard_kernel_sparse.m:
%   • the EXPONENTIAL bank  (model B) — decaying exp(-τ/τ_k) anchored at τ=0,
%     one atom per timescale, + a few raised-cosine bumps on the acausal side;
%   • the RAISED-COSINE bank (models A/C) — localized bumps tiling the lags.
% Parameters below mirror hazard_kernel_sparse.m so the atoms are identical to
% what the fit actually uses. No data needed — this only draws the bases.
% ════════════════════════════════════════════════════════════════════════

fps             = 60;
kernel_past_s   = 4.0;    % seconds of causal history   (match the main script)
kernel_future_s = 2.0;    % seconds of acausal control
nb_kernel       = 12;     % raised-cosine bumps over causal lags
nb_acausal      = 6;      % bumps over the acausal lags
bump_width      = 1.25;   % raised-cosine half-width (in centre-spacing units)
n_exp           = 12;     % exponential atoms in the bank

% ── lag grid (identical to the main script) ──────────────────────────────
back_fr = round(kernel_past_s   * fps);
fwd_fr  = round(kernel_future_s * fps);
lag_fr  = (-fwd_fr : back_fr)';
t_rel   = -lag_fr / fps;          % time rel. to un-freeze (s): past < 0 < future

% ── build the two bases exactly as hazard_kernel_sparse.m does ───────────
tau_grid_s = logspace(log10(0.05), log10(kernel_past_s), n_exp);   % 50 ms → window
[Bexp, tau_grid, is_exp] = exp_basis(lag_fr, tau_grid_s, fps, nb_acausal, bump_width);
Bflex = raised_cosine_basis(lag_fr, nb_kernel + nb_acausal, bump_width);
% raw bumps (NOT unit-normed): every bump peaks at 1, so the tiling reads cleanly.
% (the fit itself unit-norms these — see hazard_kernel_sparse.m — but that only
%  rescales edge/interior peaks; the shapes are identical.)

% ════════════════════════════════════════════════════════════════════════
% Plot
% ════════════════════════════════════════════════════════════════════════
figure('Color','w','Position',[80 80 1100 440]);

% ── (1) exponential bank ─────────────────────────────────────────────────
subplot(1,2,1); hold on
cmapE = cool(n_exp);                       % colour atoms by timescale (fast→slow)
for k = 1:n_exp
    plot(t_rel, Bexp(:,k), 'Color', cmapE(k,:), 'LineWidth', 1.6);
end
for k = n_exp+1 : size(Bexp,2)             % acausal control bumps
    plot(t_rel, Bexp(:,k), '--', 'Color', [0.6 0.6 0.6], 'LineWidth', 1);
end
xline(0, 'k--', 'LabelVerticalAlignment','bottom');
% annotate fastest / slowest timescale
text(-0.15, max(Bexp(:,1)),   sprintf('  \\tau=%.2f s', tau_grid(1)),   'Color', cmapE(1,:),   'FontSize',9);
text(-kernel_past_s*0.8, Bexp(find(lag_fr==round(kernel_past_s*fps*0.8),1), n_exp)+0.02, ...
     sprintf('\\tau=%.1f s', tau_grid(end)), 'Color', cmapE(end,:), 'FontSize',9);
xlabel('time relative to un-freeze (s)'); ylabel('basis weight')
xlim([min(t_rel) max(t_rel)]); grid on; box off

apply_generic(gca)
grid off

% ── (2) raised-cosine bank ───────────────────────────────────────────────
subplot(1,2,2); hold on
cmapR = parula(size(Bflex,2));             % colour bumps by centre position
for k = 1:size(Bflex,2)
    plot(t_rel, Bflex(:,k), 'Color', cmapR(k,:), 'LineWidth', 1.6);
end
xline(0, 'k--', 'LabelVerticalAlignment','bottom');
xlabel('time relative to un-freeze (s)'); ylabel('basis weight')
xlim([min(t_rel) max(t_rel)]); grid on; box off

apply_generic(gca)
grid off

% ════════════════════════════════════════════════════════════════════════
% Local functions (verbatim copies from hazard_kernel_sparse.m)
% ════════════════════════════════════════════════════════════════════════
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
