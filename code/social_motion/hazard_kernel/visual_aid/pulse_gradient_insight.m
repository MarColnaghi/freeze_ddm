clearvars; close all; clc;
% ════════════════════════════════════════════════════════════════════════
% pulse_gradient_insight — VISUAL AID (conceptual)
%
% One toy freeze: a strong PULSE of motion 4 s before the break, then a
% CONSTANT FLUX until it breaks. Question: where does this freeze push the
% kernel? Answer, shown here: the kernel is moved by the per-interval score
%     ∂ℓ/∂β(τ) = Σ_i (y_i − μ_i) · s_i(τ),
% i.e. every at-risk interval "votes" on β(τ) with the sign of its residual:
%   • the BREAK interval (y=1) votes UP where its history is high,
%   • every NON-BREAK interval (y=0) votes DOWN where its history is high.
% Because the window slides, the pulse sits at a DIFFERENT lag in every
% interval — and lag 4 s is the ONLY lag where it coincides with a break, not
% a survival. So the votes sum to a clean bump at 4 s.
% ════════════════════════════════════════════════════════════════════════

% ── toy parameters ───────────────────────────────────────────────────────────
fps        = 60;
dt_frames  = 6;                 % sample every 0.1 s
kernel_past_s = 4.0;  back_fr = round(kernel_past_s*fps);
entry_fr   = 30;                % 0.5 s risk-set entry
dur_s      = 5.0;   dur = round(dur_s*fps);     % freeze length; break at frame dur
pulse_u    = 4.0;               % pulse is 4 s before the break
pulse_amp  = 4.0;   pulse_wid = 0.12;           % strong, brief
flux_amp   = 1.0;               % constant flux over the last 4 s

% motion as a function of u = seconds-before-break (u ≥ 0)
motion = @(u)  flux_amp .* (u >= 0 & u <= 4) ...
             + pulse_amp .* exp(-((u - pulse_u)/pulse_wid).^2);

% ── build the at-risk intervals + their sliding-window histories ─────────────
lag_fr = (0:back_fr)';  lag_s = lag_fr/fps;  t_rel = -lag_s;   % causal window, past<0
grid   = (entry_fr:dt_frames:dur)';  if grid(end) < dur, grid=[grid;dur]; end
nI     = numel(grid);
S = zeros(nI, numel(lag_fr));  y = zeros(nI,1);  u_int = zeros(nI,1);
for i = 1:nI
    f = grid(i);  u_int(i) = (dur - f)/fps;              % this interval's time-before-break
    S(i,:) = motion(u_int(i) + lag_s');                 % history at each lag
    y(i)   = double(f == dur);                          % break on the last interval
end

% ── the gradient each interval contributes at the null (β=0, μ = base rate) ──
p0 = mean(y);                                           % base hazard; residuals sum to 0
resid = y - p0;                                         % +(1-p0) at break, -p0 elsewhere
g          = (resid' * S)';                             % net score  Σ_i (y_i-μ)s_i(τ)
g_break    = ((y==1).*resid)' * S;  g_break = g_break'; % the break interval's vote
g_nonbreak = ((y==0).*resid)' * S;  g_nonbreak = g_nonbreak';  % all survivals' votes

% ════════════════════════════════════════════════════════════════════════
% figure
% ════════════════════════════════════════════════════════════════════════
fh = figure('Color','w','Position',[60 60 900 900]);
tl = tiledlayout(3,1,'TileSpacing','compact','Padding','compact'); %#ok<NASGU>
redmap = cbrewer2('Reds', []);

% ── (1) the bout: motion vs time-before-break ────────────────────────────────
ax1 = nexttile; hold(ax1,'on');
uu = linspace(0, dur_s, 1000);
plot(ax1, -uu, motion(uu), 'k-', 'LineWidth', 2);
plot(ax1, -u_int, motion(u_int), 'o', 'MarkerSize',5, ...
    'MarkerFaceColor',[0.55 0.7 0.9], 'MarkerEdgeColor','k');   % at-risk samples
xline(ax1, 0, 'k--'); xline(ax1, -pulse_u, 'r-.', 'LineWidth',1.2);
text(ax1, -pulse_u, pulse_amp*1.02, '  pulse (−4 s)', 'Color','r', 'FontSize',13,'VerticalAlignment','bottom');
text(ax1, -2, flux_amp+0.25, 'constant flux', 'Color',[0.3 0.3 0.3], 'FontSize',13,'HorizontalAlignment','center');
text(ax1, 0, pulse_amp*0.6, ' BREAK', 'Color','k','FontWeight','bold','FontSize',13);
xlim(ax1,[-dur_s 0.25]); ylim(ax1,[-0.3 pulse_amp*1.15]);
xlabel(ax1,'time relative to break (s)'); ylabel(ax1,'social motion');
title(ax1,'one freeze:  strong pulse 4 s before the break, then constant flux','FontWeight','normal');
apply_generic(ax1,'font_size',14);

% ── (2) its sliding-window design S: the pulse traces a diagonal ─────────────
ax2 = nexttile; hold(ax2,'on');
[t_ax, ord] = sort(t_rel);
imagesc(ax2, t_ax, 1:nI, S(:,ord), 'AlphaData', ones(nI,numel(lag_fr)));
set(ax2,'YDir','reverse'); colormap(ax2, redmap); caxis(ax2,[0 pulse_amp]);
% outcome markers on the right edge
xr = max(t_ax) + 0.12;
for i=1:nI
    c = [0.9 0.9 0.87]; if y(i)==1, c=[0.94 0.63 0.15]; end
    plot(ax2, xr, i, 's', 'MarkerSize',9,'MarkerFaceColor',c,'MarkerEdgeColor',[0.5 0.5 0.5]);
end
text(ax2, xr, 0.2, 'y', 'HorizontalAlignment','center','FontWeight','bold','FontSize',13);
xlim(ax2,[min(t_ax) xr+0.15]); ylim(ax2,[0.5 nI+0.5]);
xlabel(ax2,'lag  (time relative to that interval, s)'); ylabel(ax2,'at-risk interval');
title(ax2,'sliding window: the pulse sits at a different lag in every interval (diagonal)','FontWeight','normal');
set(ax2,'YTick',[]);
apply_generic(ax2,'font_size',14); set(ax2,'YTick',[]);
cb = colorbar(ax2); cb.Label.String='motion'; cb.FontSize=12;

% ── (3) the votes sum to a bump at 4 s ───────────────────────────────────────
ax3 = nexttile; hold(ax3,'on');
plot(ax3, t_ax, g_break(ord),    '-', 'Color',[0.20 0.55 0.25],'LineWidth',1.6);
plot(ax3, t_ax, g_nonbreak(ord), '-', 'Color',[0.80 0.30 0.25],'LineWidth',1.6);
plot(ax3, t_ax, g(ord),          'k-','LineWidth',2.6);
plot(ax3, t_ax, 0*t_ax, 'k-','LineWidth',0.6);
[~,ip] = max(g); xline(ax3, t_rel(ip), 'r-.','LineWidth',1.2);
xlim(ax3,[min(t_ax) 0]);
xlabel(ax3,'lag  \tau  (time relative to break, s)'); ylabel(ax3,'push on  \beta(\tau)');
legend(ax3, {'break interval votes (y=1, up)','survival intervals vote (y=0, down)','NET push'}, ...
    'Location','northwest','Box','off','FontSize',12);
apply_generic(ax3,'font_size',14);
title(ax3,'net push = sum of (residual x history) over intervals  -  peaks at the pulse-to-break lag (4 s)', ...
    'FontWeight','normal');

set([ax1 ax2 ax3],'Layer','top');
ax1.Toolbar.Visible='off'; ax2.Toolbar.Visible='off'; ax3.Toolbar.Visible='off';
outdir = fileparts(mfilename('fullpath'));
exportgraphics(fh, fullfile(outdir,'pulse_gradient_insight.png'), 'Resolution',200);   % PNG preview
paths = path_generator('folder','social_motion/hazard_kernel/visual_aid','imfirst',false);
exporter(fh, paths, 'pulse_gradient_insight.pdf');   % vector PDF -> paths.fig (repo convention)
fprintf('peak of net push at t_rel = %.2f s\n', t_rel(ip));
fprintf('saved: %s\n       %s\n', fullfile(outdir,'pulse_gradient_insight.png'), fullfile(paths.fig,'pulse_gradient_insight.pdf'));
