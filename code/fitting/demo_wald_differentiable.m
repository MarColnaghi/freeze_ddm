%% Why the Wald likelihood is DIFFERENTIABLE (and the simulator's is not)
%   Constant-drift, no-leak DDM. The first-passage density is the Wald
%   (inverse Gaussian). We differentiate the log-likelihood in mu three ways
%   and show they agree:
%     (1) the closed-form score   dL/dmu = sum (a - mu t)/sigma^2   [by hand]
%     (2) finite differences of the analytic L(mu)                  [numeric check]
%   ...and contrast with the MONTE-CARLO estimate of the same L(mu), which is
%   a jagged staircase with no usable gradient.
%
%   The point: "differentiable" is a property of the LIKELIHOOD (the density,
%   averaged over noise), not of the model. A single path crosses the bound at
%   a discrete instant -- a hard threshold -- but averaging that indicator
%   against the continuous Gaussian noise smooths it into a C^inf function of
%   the parameters. The Wald IS that smoothed object in closed form.
%
%   Run headless:
%     matlab -batch "run('fitting/demo_wald_differentiable.m')"
clearvars; clc;
this_dir = fileparts(mfilename('fullpath'));
code_dir = fileparts(this_dir);
addpath(fullfile(code_dir,'sims','simulators'));

sigma = 1.0;  a = 1.2;  mu0 = 1.0;      % truth
rng(3);

% ---- a fixed dataset: exact Wald draws at the true mu -----------------------
Ndat = 500;
m_ig = a/mu0;  lam_ig = a^2/sigma^2;    % inverse-Gaussian mean & shape
nu = randn(Ndat,1); y = nu.^2;
rt = m_ig + (m_ig^2*y)/(2*lam_ig) - (m_ig/(2*lam_ig))*sqrt(4*m_ig*lam_ig*y + m_ig^2*y.^2);
u  = rand(Ndat,1);
sw = u > m_ig./(m_ig+rt);  rt(sw) = m_ig^2 ./ rt(sw);   % Michael-Schucany-Haas

% ---- (A) the density is a smooth family in mu -------------------------------
mus_show = [0.6 0.9 1.2 1.5];
tt = linspace(0.01, 5, 800)';
Fd = zeros(numel(tt), numel(mus_show));
for j = 1:numel(mus_show), Fd(:,j) = wald_fpt(tt, a, mus_show(j), sigma); end

% ---- (B) analytic log-likelihood L(mu) and its closed-form gradient ---------
mgrid = linspace(0.45, 1.65, 241)';
ellA  = arrayfun(@(mm) sum(log(max(wald_fpt(rt,a,mm,sigma),realmin))), mgrid);
scoreA = arrayfun(@(mm) sum((a - mm*rt)/sigma^2), mgrid);   % closed form dL/dmu
dell_fd = gradient(ellA, mgrid);                            % numeric derivative

% ---- (C) Monte-Carlo estimate of the SAME L(mu): jagged, no gradient --------
%   Estimate f(.|mu) by simulating paths and smoothing (ksdensity), then score
%   the fixed data. INDEPENDENT seeds across the grid -> the finite-sample noise
%   is visible: the estimate jitters around the analytic curve and has no clean
%   slope you could differentiate.
dt = 1/120;  Nmc = 3000;  Ts = round(6/dt);
mgrid_mc = linspace(0.55, 1.5, 49)';
ellMC = nan(size(mgrid_mc));
for i = 1:numel(mgrid_mc)
    mm = mgrid_mc(i);
    drift = mm*ones(Ts,1);
    rk = nan(Nmc,1);
    for s = 1:Nmc
        k = sim_leaky_accumulator(drift, a, sigma, 0, dt, (i-1)*Nmc + s, 0);
        if ~isnan(k), rk(s) = k*dt; end
    end
    good = rk(~isnan(rk));
    fpdf = ksdensity(good, rt);          % smooth MC density at the data points
    ellMC(i) = sum(log(max(fpdf, realmin)));
end

in = 2:numel(mgrid)-1;   % interior: gradient() is one-sided at the two endpoints
fprintf('max |closed-form score - finite-diff of L| (interior) = %.3e\n', ...
        max(abs(scoreA(in) - dell_fd(in))));
[~,imax] = max(ellA);
fprintf('analytic MLE mu-hat = %.4f  (truth %.2f)\n', mgrid(imax), mu0);

%% ---- figure ---------------------------------------------------------------
figure('Color','w','Position',[40 40 1550 470]);
tiledlayout(1,3,'TileSpacing','compact','Padding','compact');
cmap = [0.15 0.35 0.65; 0.20 0.60 0.55; 0.85 0.55 0.15; 0.75 0.25 0.25];

% (A)
ax = nexttile; hold(ax,'on')
for j = 1:numel(mus_show)
    plot(ax, tt, Fd(:,j), '-', 'LineWidth', 2.2, 'Color', cmap(j,:), ...
         'DisplayName', sprintf('\\mu = %.1f', mus_show(j)));
end
xlabel(ax,'first-passage time t (s)'); ylabel(ax,'f(t \mid \mu)');
title(ax,'A. the density slides SMOOTHLY with \mu','FontWeight','normal');
legend(ax,'Location','northeast','Box','off'); xlim(ax,[0 5]); set(ax,'FontSize',13)
box(ax,'off')

% (B)
ax = nexttile; hold(ax,'on')
plot(ax, mgrid, ellA, '-', 'LineWidth', 2.6, 'Color', [0.15 0.35 0.65], ...
     'DisplayName','analytic  L(\mu) = \Sigma log f(t_i\mid\mu)');
plot(ax, mgrid_mc, ellMC, 'o-', 'LineWidth', 1.0, 'MarkerSize', 3.5, ...
     'Color', [0.85 0.45 0.15], 'MarkerFaceColor',[0.85 0.45 0.15], ...
     'DisplayName','Monte-Carlo estimate (N=3000)');
xline(ax, mu0, 'k:', 'truth', 'LineWidth', 1.2, 'FontSize', 11, 'HandleVisibility','off');
xlabel(ax,'\mu'); ylabel(ax,'log-likelihood');
title(ax,'B. analytic L is smooth; MC L is jagged','FontWeight','normal');
legend(ax,'Location','south','Box','off','FontSize',10);
set(ax,'FontSize',13); box(ax,'off')
lo = min([ellA; ellMC(isfinite(ellMC))]);
ylim(ax,[lo-5, max(ellA)+8])

% (C)
ax = nexttile; hold(ax,'on')
plot(ax, mgrid, scoreA, '-', 'LineWidth', 3.0, 'Color', [0.75 0.80 0.90], ...
     'DisplayName','closed form  \Sigma(a-\mu t_i)/\sigma^2');
plot(ax, mgrid(1:6:end), dell_fd(1:6:end), 'o', 'MarkerSize', 5, ...
     'Color',[0.15 0.35 0.65], 'MarkerFaceColor',[0.15 0.35 0.65], ...
     'DisplayName','finite-diff of L(\mu)');
yline(ax, 0, 'k:', 'LineWidth', 1, 'HandleVisibility','off');
xline(ax, mgrid(imax), '-', 'MLE', 'Color',[0.4 0.4 0.4], 'LineWidth',1, 'FontSize',11, 'HandleVisibility','off');
xlabel(ax,'\mu'); ylabel(ax,'dL/d\mu');
title(ax,'C. the gradient is a formula (\equiv finite diff)','FontWeight','normal');
legend(ax,'Location','northeast','Box','off','FontSize',10);
set(ax,'FontSize',13); box(ax,'off')

exportgraphics(gcf, fullfile(this_dir,'demo_wald_differentiable.png'), 'Resolution', 140);
fprintf('done.\n');
