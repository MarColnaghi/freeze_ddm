clearvars; close all; clc;
% ════════════════════════════════════════════════════════════════════════
% DEMO: the BOUND (not the leak) is what bends the recovered kernel — and
% you can't get rid of the bend by re-aligning; you can only MOVE it.
%
% All three panels use the project's OWN first-passage engine,
% sim_leaky_accumulator, with lambda = 0 (PERFECT accumulator, NO leak). A
% perfect accumulator weights past evidence EQUALLY — its internal kernel is
% FLAT. What changes across panels is the bound and the alignment:
%
%   Panel A — bounded, ONSET-aligned    -> PRIMACY  (declines from onset)
%   Panel B — bounded, RESPONSE-aligned -> RECENCY  (rises to the un-freeze)
%   Panel C — UNBOUNDED, fixed horizon  -> FLAT
%
% Same perfect integrator throughout. A & B are the SAME bounded bouts, only
% the alignment differs: the absorbing bound bends the kernel at BOTH ends
% (primacy from onset, recency to the response). Only REMOVING the bound
% (Panel C) gives a flat kernel — which is the true internal weighting.
% ════════════════════════════════════════════════════════════════════════

here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..', '..', '..', 'sims', 'simulators'));   % sim_leaky_accumulator
rng(7);

% ── shared knobs ─────────────────────────────────────────────────────────
fps    = 60;   dt = 1/fps;
L      = 100;    % past lags shown in the kernel (frames)
sigma  = 1.0;    % diffusion SD inside the accumulator (Euler–Maruyama)
mu0    = 2.2;    % baseline drift RATE (urge to un-freeze) -- tuned so bouts cross
                 %   AROUND the window L, so A sees attrition and B still has events
beta   = 4.2;    % gain on the social-motion fluctuations we reverse-correlate
lam    = 5;      % small ridge to denoise the raw per-lag fit (shape unaffected)
LAMBDA = 0;      % LEAK RATE = 0  -> PERFECT integrator (flat internal kernel)

% ════════════════════════════════════════════════════════════════════════
% BOUNDED perfect integrator (Panels A & B share these bouts)
% ════════════════════════════════════════════════════════════════════════
nBout  = 50000;  % freeze bouts
theta  = 3.5;    % un-freeze threshold (the BOUND)
maxLen = 300;    % integration cap (frames)

Xc = cell(nBout,1); Yc = cell(nBout,1);
Mo     = zeros(nBout, L);                 % ONSET-aligned motion (frames 1..L from onset)
kcross = nan(nBout, 1);                   % un-freeze frame per bout (NaN = never crossed)
for b = 1:nBout
    m     = randn(maxLen,1);              % social-motion fluctuations (z-scored)
    drive = mu0 + beta*m;                 % per-frame drift RATE
    k = sim_leaky_accumulator(drive, theta, sigma, LAMBDA, dt, 1000+b);   % first-passage

    Mo(b,:)   = m(1:L).';                 % onset-aligned predictors (ALL bouts; frame 1 = onset)
    kcross(b) = k;                        % un-freeze frame (for the onset-aligned outcome)

    if isnan(k) || k <= L, continue; end  % response-aligned needs full L-lag history
    tt = (L+1:k).';                       % at-risk steps with complete history (nstep x 1)
    H  = zeros(numel(tt), L);             % history matrix, guaranteed nstep x L
    for j = 1:L
        H(:,j) = m(tt - j);               % column j = motion j frames before the step
    end
    Xc{b} = H;                            % lag 1 = most recent
    Yc{b} = double(tt == k);              % y=1 only on the crossing step
end
keep = ~cellfun(@isempty, Xc);
X1 = vertcat(Xc{keep});
y1 = vertcat(Yc{keep});

% ── regime check: does the crossing distribution STRADDLE the window L? ───
% Need mass BELOW L (attrition -> Panel-A primacy) AND above L (Panel-B events).
kf = kcross(~isnan(kcross));
fprintf(['Crossing times: median=%.0f fr | within window (<=L=%d): %.0f%% ', ...
         '| beyond L: %.0f%% | never crossed: %.0f%%\n'], ...
        median(kf), L, 100*mean(kf<=L), 100*mean(kf>L), 100*mean(isnan(kcross)));

% ── (B) RESPONSE-aligned hazard fit — expect RECENCY ─────────────────────
b1     = logit_fit(X1, y1, lam);          % hazard GLM (logistic on person-period rows)
kern_b = b1(2:end);                       % kernel over lags 1..L BEFORE the un-freeze
fprintf('Panel B response-aligned (bounded, theta=%.2f): %d rows, %d events\n', ...
        theta, numel(y1), sum(y1));

% ── (A) ONSET-aligned fit on the SAME bouts — expect PRIMACY ─────────────
% Outcome = "un-freeze within a horizon"; predictors = motion frames 1..L from
% ONSET. Onset is exogenous (NOT selected by the bound), so this is not locked
% to the decision — yet the absorbing bound still bends the kernel, the OTHER
% way: late-from-onset frames act only on bouts still frozen, so their leverage
% fades => primacy (declining from onset), NOT flat.
horizon = round(median(kcross(~isnan(kcross))));   % balances the outcome
y_on    = double(kcross <= horizon);               % escaped by horizon (NaN -> 0)
b_on    = logit_fit(Mo, y_on, lam);
kern_on = b_on(2:end);                             % kernel over frames 1..L FROM onset
fprintf('Panel A onset-aligned   (bounded, horizon=%d fr): %d bouts, %.0f%% escaped-by-horizon\n', ...
        horizon, nBout, 100*mean(y_on));

% ════════════════════════════════════════════════════════════════════════
% (C) UNBOUNDED, fixed horizon — expect FLAT
% ════════════════════════════════════════════════════════════════════════
nTrial = 40000;  % trials
Thor   = L;      % fixed horizon = kernel window
a      = 0.9;    % choice sensitivity to the accumulated evidence

Mtr = randn(nTrial, Thor);                % each row: Thor frames of motion
v   = nan(nTrial,1);
for i = 1:nTrial
    drive = mu0 + beta*Mtr(i,:)';
    [~, traj] = sim_leaky_accumulator(drive, Inf, sigma, LAMBDA, dt, 500000+i); % never stops
    v(i) = traj(end);                     % accumulator state at the fixed horizon
end
p  = 1 ./ (1 + exp(-a*(v - median(v))));  % ONE decision, no absorption
y2 = double(rand(nTrial,1) < p);

b2     = logit_fit(Mtr, y2, lam);         % reverse-correlation: choice on each frame
kern_f = flipud(b2(2:end));               % re-index to "lag from decision" (1=most recent)
fprintf('Panel C unbounded       (theta=Inf): %d trials, %.0f%% "break"\n', ...
        nTrial, 100*mean(y2));

% ── RAW recovered kernels (log-odds units, no normalisation) ─────────────
lag  = (1:L)';                            % response/decision lag: 1 = most recent
fro  = (1:L)';                            % frames from onset:      1 = onset

s_on = polyfit(fro, kern_on, 1);
s_b  = polyfit(lag, kern_b,  1);
s_f  = polyfit(lag, kern_f,  1);
fprintf('\nKernel slope  [0 = flat] :\n');
fprintf('   A onset-aligned  vs frames-from-onset : % .5f   <-- primacy (declines from onset)\n', s_on(1));
fprintf('   B response-align vs lag-before-event  : % .5f   <-- recency (rises to the event)\n',   s_b(1));
fprintf('   C unbounded      vs lag               : % .5f   <-- flat\n',                            s_f(1));

% ════════════════════════════════════════════════════════════════════════
% FIGURE — 3 panels
% ════════════════════════════════════════════════════════════════════════
figure('Color','w','Position',[60 80 1500 420]);

% (A) onset-aligned, bounded -> primacy
subplot(1,3,1); hold on
plot(fro, kern_on, '-o', 'Color',[0.55 0.25 0.6], 'LineWidth',1.8, 'MarkerSize',3, ...
     'MarkerFaceColor',[0.55 0.25 0.6]);
yline(mean(kern_on), 'k--', 'LineWidth',1);
xlabel('time from freeze ONSET (frames)   1 = onset');
ylabel('recovered kernel weight (log-odds)');
title({'\bfA  —  bounded, ONSET-aligned','\rm\theta finite, \lambda=0  \Rightarrow  PRIMACY'});
xlim([1 L]); box off; grid on

% (B) response-aligned, bounded -> recency
subplot(1,3,2); hold on
plot(lag, kern_b, '-o', 'Color',[0.85 0.2 0.2], 'LineWidth',1.8, 'MarkerSize',3, ...
     'MarkerFaceColor',[0.85 0.2 0.2]);
yline(mean(kern_b), 'k--', 'LineWidth',1);
xlabel('lag before UN-FREEZE (frames)   1 = most recent');
ylabel('recovered kernel weight (log-odds)');
title({'\bfB  —  bounded, RESPONSE-aligned','\rm\theta finite, \lambda=0  \Rightarrow  RECENCY'});
xlim([1 L]); box off; grid on

% (C) unbounded -> flat
subplot(1,3,3); hold on
plot(lag, kern_f, '-o', 'Color',[0.2 0.4 0.85], 'LineWidth',1.8, 'MarkerSize',3, ...
     'MarkerFaceColor',[0.2 0.4 0.85]);
yline(mean(kern_f), 'k--', 'LineWidth',1);
xlabel('lag from decision (frames)   1 = most recent');
ylabel('recovered kernel weight (log-odds)');
title({'\bfC  —  UNBOUNDED, fixed horizon','\rm\theta=\infty, \lambda=0  \Rightarrow  FLAT'});
xlim([1 L]); box off; grid on

sgtitle('Same perfect integrator (\lambda=0): the BOUND bends the kernel at BOTH ends — only removing it gives FLAT');

% ════════════════════════════════════════════════════════════════════════
% Local function: ridge-penalised logistic regression via IRLS
% ════════════════════════════════════════════════════════════════════════
function b = logit_fit(X, y, lam)
% Logistic regression with a small L2 ridge (intercept unpenalised).
% Returns b = [intercept; per-column weights].
    [n, p] = size(X);
    X = [ones(n,1), X];
    b = zeros(p+1, 1);
    R = lam * eye(p+1); R(1,1) = 0;           % no penalty on the intercept
    for it = 1:100
        eta = min(max(X*b, -30), 30);
        mu  = 1 ./ (1 + exp(-eta));
        w   = max(mu .* (1 - mu), 1e-6);
        z   = eta + (y - mu) ./ w;            % IRLS working response
        bn  = (X' * (w .* X) + R) \ (X' * (w .* z));
        if max(abs(bn - b)) < 1e-8, b = bn; break; end
        b = bn;
    end
end
