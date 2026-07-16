clearvars; close all; clc;
% ════════════════════════════════════════════════════════════════════════
% irls_iterations — VISUAL AID (conceptual)
%
% How the logistic (hazard) GLM is fit: iteratively reweighted least squares
% = Newton's method on the concave Bernoulli log-likelihood ℓ(b). At the
% current estimate b^(t) the method builds a QUADRATIC SURROGATE
%     q_t(b) = ℓ(b^(t)) + ℓ'(b^(t))(b−b^(t)) + ½ ℓ''(b^(t))(b−b^(t))²,
% which is exactly the weighted least-squares objective (fit4: the update
% b ← (XᵀWX)⁻¹XᵀWz maximises this parabola). The next iterate b^(t+1) is the
% PEAK of that parabola. Repeating climbs to the MLE b*.
%
% Two panels:
%   (L) b-space  : ℓ(b) with two surrogates q_t, q_{t+1} and the jumps
%                  b^(t) → b^(t+1) → b^(t+2) → b*.
%   (R) data-space: the fitted hazard σ(b·x) sharpening toward the data.
% ════════════════════════════════════════════════════════════════════════

rng(3);

% ── toy 1-parameter logistic:  h = σ(b·x) ────────────────────────────────
n     = 60;
x     = linspace(-3, 3, n)';
btrue = 1.3;
y     = double(rand(n,1) < 1./(1+exp(-btrue*x)));

sig  = @(t) 1 ./ (1 + exp(-t));
LL   = @(b) sum( y.*log(sig(b*x)+eps) + (1-y).*log(1-sig(b*x)+eps) );  % ℓ(b)
grad = @(b) sum( (y - sig(b*x)).*x );                                  % ℓ'(b)
hess = @(b) -sum( sig(b*x).*(1-sig(b*x)).*x.^2 );                      % ℓ''(b) < 0

% ── MLE b* (many Newton steps) and two iterations from a start b0 ─────────
bstar = 0; for it = 1:100, bstar = bstar - grad(bstar)/hess(bstar); end
b0 = 0.20;
b1 = b0 - grad(b0)/hess(b0);
b2 = b1 - grad(b1)/hess(b1);
fprintf('b0=%.3f  b1=%.3f  b2=%.3f  b*=%.3f\n', b0, b1, b2, bstar);

% surrogate parabola at expansion point bt (vertex sits at next iterate)
q = @(b,bt) LL(bt) + grad(bt).*(b-bt) + 0.5*hess(bt).*(b-bt).^2;

cA = [0.20 0.45 0.80];   % surrogate q_t   (blue)
cB = [0.85 0.45 0.10];   % surrogate q_{t+1} (orange)
cS = [0.15 0.60 0.25];   % MLE (green)

% ════════════════════════════════════════════════════════════════════════
fh = figure('Color','w','Position',[60 60 1150 480]);
tl = tiledlayout(1,2,'TileSpacing','compact','Padding','compact'); %#ok<NASGU>

% ── (L) b-space: log-likelihood + two quadratic surrogates ───────────────
axA = nexttile; hold(axA,'on');
bb  = linspace(min([b0 bstar])-0.5, max([b0 bstar])+1.0, 500)';
ll  = arrayfun(LL, bb);
plot(axA, bb, ll, 'k-', 'LineWidth', 2.4);                    % true ℓ(b)
plot(axA, bb, q(bb,b0), '--', 'Color', cA, 'LineWidth', 1.8); % surrogate q_t
plot(axA, bb, q(bb,b1), '--', 'Color', cB, 'LineWidth', 1.8); % surrogate q_{t+1}
ylo = min(ll) - 0.06*range(ll); yhi = max(ll) + 0.10*range(ll);
ylim(axA,[ylo yhi]); xlim(axA,[bb(1) bb(end)]);

% points on the true curve at each iterate
plot(axA, b0, LL(b0), 'ks', 'MarkerFaceColor','k','MarkerSize',9);
plot(axA, b1, LL(b1), 'ko', 'MarkerFaceColor','w','MarkerSize',8,'LineWidth',1.5);
plot(axA, b2, LL(b2), 'ko', 'MarkerFaceColor','w','MarkerSize',8,'LineWidth',1.5);
% surrogate peaks (= next iterate)
plot(axA, b1, q(b1,b0), 'o', 'Color',cA,'MarkerFaceColor',cA,'MarkerSize',7);
plot(axA, b2, q(b2,b1), 'o', 'Color',cB,'MarkerFaceColor',cB,'MarkerSize',7);
% vertical guides down to the axis showing the new iterate
plot(axA,[b1 b1],[ylo q(b1,b0)],':','Color',cA,'LineWidth',1.2);
plot(axA,[b2 b2],[ylo q(b2,b1)],':','Color',cB,'LineWidth',1.2);
xline(axA, bstar, '-', 'Color', cS, 'LineWidth', 1.6);

% Newton jumps as arrows along a low y-level
ylev = ylo + 0.035*range(ll);
quiver(axA, b0, ylev, b1-b0, 0, 0, 'Color',cA,'LineWidth',1.8,'MaxHeadSize',0.6);
quiver(axA, b1, ylev, b2-b1, 0, 0, 'Color',cB,'LineWidth',1.8,'MaxHeadSize',0.6);

text(axA, b0, ylo-0.02*range(ll), 'b^{(t)}',   'Color','k', 'HorizontalAlignment','center','VerticalAlignment','top','FontSize',13);
text(axA, b1, ylo-0.02*range(ll), 'b^{(t+1)}', 'Color',cA,'HorizontalAlignment','center','VerticalAlignment','top','FontSize',13);
text(axA, b2, ylo-0.02*range(ll), 'b^{(t+2)}', 'Color',cB,'HorizontalAlignment','center','VerticalAlignment','top','FontSize',13);
text(axA, bstar, yhi, ' b*', 'Color',cS,'HorizontalAlignment','left','VerticalAlignment','top','FontSize',13,'FontWeight','bold');
text(axA, b0+0.03, q(b0+0.03,b0)-0.01*range(ll), ' q_t', 'Color',cA,'FontSize',12);

xlabel(axA,'coefficient  b'); ylabel(axA,'log-likelihood  \ell(b)');
title(axA,'each step maximises a quadratic (weighted-LS) surrogate','FontWeight','normal');
set(axA,'Box','off','FontSize',13,'Layer','top','TickDir','out');

% ── (R) data-space: fitted hazard sharpening toward the data ─────────────
axB = nexttile; hold(axB,'on');
xx = linspace(-3.2, 3.2, 300)';
plot(axB, x(y==0), zeros(sum(y==0),1)+0.02, 'o', 'MarkerSize',5,'MarkerEdgeColor',[0.5 0.5 0.5],'MarkerFaceColor','w');
plot(axB, x(y==1), ones(sum(y==1),1)-0.02,  'o', 'MarkerSize',5,'MarkerEdgeColor',[0.5 0.5 0.5],'MarkerFaceColor',[0.8 0.8 0.8]);
plot(axB, xx, sig(b0*xx), '-', 'Color',[0.72 0.72 0.72],'LineWidth',1.8);
plot(axB, xx, sig(b1*xx), '--','Color',cA,'LineWidth',1.9);
plot(axB, xx, sig(b2*xx), '--','Color',cB,'LineWidth',1.9);
plot(axB, xx, sig(bstar*xx),'-','Color',cS,'LineWidth',2.4);
ylim(axB,[-0.06 1.06]); xlim(axB,[-3.2 3.2]);
xlabel(axB,'feature  x_i'); ylabel(axB,'predicted hazard  h = \sigma(b\cdotx_i)');
title(axB,'the fitted hazard sharpens toward the data each step','FontWeight','normal');
legend(axB, {'y = 0','y = 1','b^{(t)}','b^{(t+1)}','b^{(t+2)}','b* (converged)'}, ...
    'Location','southeast','Box','off','FontSize',11);
set(axB,'Box','off','FontSize',13,'Layer','top','TickDir','out');

% ── export ───────────────────────────────────────────────────────────────
outdir = fileparts(mfilename('fullpath'));
exportgraphics(fh, fullfile(outdir,'irls_iterations.png'), 'Resolution',200);   % PNG preview
paths = path_generator('folder','social_motion/hazard_kernel/visual_aid','imfirst',false);
exporter(fh, paths, 'irls_iterations.pdf');   % vector PDF -> paths.fig (repo convention)
fprintf('saved: %s\n       %s\n', fullfile(outdir,'irls_iterations.png'), fullfile(paths.fig,'irls_iterations.pdf'));
