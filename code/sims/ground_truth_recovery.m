clearvars; close all; clc;
% ════════════════════════════════════════════════════════════════════════
% GROUND-TRUTH RECOVERY — single leaky-accumulator continuum.
%
% One model (sim_leaky_accumulator), swept over the integration time-constant τ:
%     τ = ∞  (λ=0)        → perfect accumulator   (flat kernel, infinite memory)
%     τ large→small (λ↑)  → leaky integrator       (exponential kernel, memory τ)
%     τ → dt  (λ→1/dt)    → instantaneous/extrema  (delta kernel, no memory)
% Driven by the REAL sm traces (same flies/onsets). Each τ's bound is calibrated
% to the real median freeze duration, and we VERIFY the truncated duration
% distributions match (KS vs real) before trusting the kernel comparison.
%
% Two readouts answer the two questions:
%   (1) TIME HORIZON  → recovered kernel timescale (τ_com) vs true τ.
%   (2) TEMPORAL-STRUCTURE SENSITIVITY → the time-shuffle gap (real − shuffled CV
%       log-lik) vs true τ: short-τ (recency) models should be hurt MORE by
%       destroying temporal order than long-τ (integrator) models.
% ════════════════════════════════════════════════════════════════════════

addpath('/Users/marcocolnaghi/PhD/freeze_ddm/code/sims/simulators');
addpath('/Users/marcocolnaghi/PhD/freeze_ddm/code/social_motion/proximity');

% ── Regime (rebalanced: baseline sets the timescale, sm is a MODULATION) ──
fps        = 60;
P.dt       = 1/fps;
P.mu0      = 10;     % baseline drive — sets the ~real-median timescale
P.beta     = 10.0;     % social gain per 1 SD of sm (moderate: modulates, doesn't dominate)
P.sigma    = 1.0;     % diffusion SD
P.seed     = 1000;
maxT_fr    = 630;     % 10.5 s integration cap
trunc_point= 30;      % 0.5 s min duration (mirror real pipeline)
n_rep      = 3;       % replications per onset (× events; raise for finer kernels)

% τ continuum: Inf (accumulator) → 0.02 s (≈ instantaneous; τ=dt=1/60 is the floor)
tau_list   = [Inf] %0.5 0.25 0.1 0.05 0.02];

n_shuffle_recovery = 20;   % time-shuffle reps for the temporal-structure gap (0 = skip)

opts = struct('fps',fps,'kernel_past_s', 6,'kernel_future_s', 2,'dt_frames', 3, ...
    'nb_kernel', 18,'nb_acausal', 6,'bump_width', 1.5,'nb_baseline',6,'cv_folds',2, ...
    'n_shuffle',n_shuffle_recovery,'shuffle_mode','full','trunc_point',trunc_point, ...
    'do_plot',true,'do_stage2',false,'use_penalty',true,'anchor_zero',true);   % penalised kernel (1-SE λ) + ACF/ETA control + τ=0 knot

% ── Load real bouts + sm, build template (real onsets/flies/conditions) ───
threshold_imm=2; threshold_mob=2; threshold_pc=4;
id_code = sprintf('imm%d_mob%d_pc%d', threshold_imm, threshold_mob, threshold_pc);
paths = path_generator('folder','social_motion/proximity','bouts_id',id_code,'imfirst',false);
bouts = importdata(fullfile(paths.dataset,'bouts.mat'));
motion_cache = importdata(fullfile(paths.cache_path,'motion_cache.mat'));
thresholds = define_thresholds('le_window', struct('le_window_sl',[5 55],'le_window_fl',[5 55]));
bouts = bouts_formatting(bouts, thresholds);
bl_t  = data_parser_new(bouts,'type','immobility','period','loom','window','le', ...
        'nloom',2:20,'min_dur',trunc_point);

if ismember('durations', bl_t.Properties.VariableNames)
    real_dur_all = double(bl_t.durations);
else
    real_dur_all = double(bl_t.ends - bl_t.onsets);
end
real_med  = median(real_dur_all);
real_durs = real_dur_all(real_dur_all > trunc_point & real_dur_all < maxT_fr);   % for KS
fprintf('Template: %d real freezes, real median = %.0f fr (%.2f s)\n', height(bl_t), real_med, real_med/fps);

% per-bout sm input (onset forward, capped), z-scored globally; NaN→baseline
nB = height(bl_t); sm_in = cell(nB,1);
for b = 1:nB
    sm = motion_cache(bl_t.fly(b)); sm = sm(:); Lf = numel(sm); on = bl_t.onsets(b);
    cap = min(maxT_fr, Lf - on + 1);
    sm_in{b} = sm(on : on + cap - 1);
end
% Scale by SD only — KEEP sm POSITIVE; do NOT center. Social motion is a
% non-negative drive that ADDS to baseline (more motion → faster un-freeze);
% absence of motion = baseline (mu0). Centering would let below-average sm push
% the accumulator backward, which is unphysical for the generative model. (The
% analysis z-scores its regression predictor, which is fine — centering a GLM
% predictor changes only the intercept, not the kernel shape.)
allsm = cell2mat(cellfun(@(v) v(:), sm_in, 'UniformOutput', false));
sm_sd = std(allsm,'omitnan');
for b = 1:nB
    sm_in{b} = sm_in{b} / sm_sd;          % scale only, stays >= 0
    sm_in{b}(isnan(sm_in{b})) = 0;        % missing sm → no social drive (baseline)
end
cap_b = cellfun(@numel, sm_in);
fprintf('sm drive scaled by SD=%.3f (kept positive, not centered)\n', sm_sd);

% ── Sweep over τ: calibrate θ → simulate → analyse ────────────────────────
R = struct('tau',{},'theta',{},'res',{},'ks',{},'med',{},'cens',{},'filt',{}, ...
           'gap',{},'pval',{},'simdur',{});
fprintf('\nCalibrating θ to median %.0f and analysing each τ ...\n', real_med);
for m = 1:numel(tau_list)
    tau = tau_list(m);
    lam = 0; if ~isinf(tau), lam = 1/tau; end

    theta = calibrate_theta(@(th) sim_leaky_all(th, lam, sm_in, P), real_med, trunc_point);

    DUR = []; CENS = [];
    for rep = 1:n_rep
        Prep = P; Prep.seed = P.seed + rep*100000;
        [d, c] = sim_leaky_all(theta, lam, sm_in, Prep);
        d(c) = cap_b(c);
        DUR = [DUR; d]; CENS = [CENS; c]; %#ok<AGROW>
    end

    % distribution-match diagnostic (truncated, vs real)
    sd = DUR(~CENS); sd = sd(sd > trunc_point & sd < maxT_fr);
    ks = NaN; if numel(sd) > 50, [~,~,ks] = kstest2(sd, real_durs); end
    filt = mean(DUR <= trunc_point);          % fraction lost to the 0.5 s filter

    bl = table();
    bl.fly          = repmat(bl_t.fly,          n_rep, 1);
    bl.onsets       = repmat(bl_t.onsets,       n_rep, 1);
    bl.moving_flies = repmat(bl_t.moving_flies, n_rep, 1);
    bl.sloom        = repmat(bl_t.sloom,        n_rep, 1);
    bl.dur_frames   = DUR; bl.is_censored = CENS;
    bl.bout_id      = (1:height(bl))';
    bl = bl(bl.dur_frames > trunc_point, :);

    if isinf(tau), opts.fig_title='accumulator (τ=∞)'; else, opts.fig_title=sprintf('τ=%.2g s',tau); end
    res = run_hazard_kernel(bl, motion_cache, opts);

    gap = NaN; pval = NaN;
    if isfield(res,'ll_real'), gap = res.ll_real - mean(res.ll_null); pval = res.pval; end

    R(m).tau=tau; R(m).theta=theta; R(m).res=res; R(m).ks=ks;
    R(m).med=median(DUR(~CENS)); R(m).cens=mean(CENS); R(m).filt=filt;
    R(m).gap=gap; R(m).pval=pval; R(m).simdur=sd;
    fprintf(['τ=%6.3g | θ=%5.2f | med=%4.0f cens=%4.1f%% filt=%4.1f%% KS=%.3f | ' ...
             'τ_com=%.2f τ_exp=%.2f | shuffleΔ=%.4f (p=%.3f)\n'], ...
        tau, theta, R(m).med, 100*R(m).cens, 100*filt, ks, ...
        res.tau_com, res.tau_hat, gap, pval);
end

taus = [R.tau]; finite = ~isinf(taus);

% ── Figure 1: recovered kernels across the continuum (absolute amplitude) ─
cmap = parula(numel(R));
figure('Color','w','Position',[80 80 820 480]); hold on
hh = gobjects(1,numel(R)); lab = cell(1,numel(R));
for m = 1:numel(R)
    res=R(m).res; cm=res.caus_mask; tr=res.t_rel(cm); kc=res.kernel(cm);
    hh(m)=plot(tr,kc,'Color',cmap(m,:),'LineWidth',1.8);
    if isinf(R(m).tau), lab{m}='τ=∞ (accum)'; else, lab{m}=sprintf('τ=%.2g',R(m).tau); end
end
yline(0,'k:'); xline(0,'--k','escape');
xlabel('Time relative to un-freeze (s)   (causal: τ ≤ 0)')
ylabel('Recovered kernel \beta(\tau)  (per SD of sm)')
title('Recovered kernels across the leaky-accumulator continuum')
legend(hh, lab, 'Location','northwest','box','off')

% ── Figure 2: TIME-HORIZON recovery — recovered τ_com/τ_exp vs true τ ─────
tcom = arrayfun(@(k) R(k).res.tau_com, 1:numel(R));
texp = arrayfun(@(k) R(k).res.tau_hat, 1:numel(R));
figure('Color','w','Position',[80 80 560 470]); hold on
plot([0.02 5],[0.02 5],'k:','LineWidth',1)                          % identity
plot(taus(finite), tcom(finite), 'o-','LineWidth',1.6,'DisplayName','\tau_{com}')
plot(taus(finite), texp(finite), 's--','LineWidth',1.2,'DisplayName','\tau_{exp}')
if any(~finite)   % accumulator as a ceiling
    yline(tcom(~finite), ':', 'accum \tau_{com}', 'Color',[0.4 0.4 0.4]);
end
set(gca,'XScale','log','YScale','log'); grid on
xlabel('True leak \tau (s)'); ylabel('Recovered \tau (s)')
title('Time-horizon recovery (identity = perfect recovery)')
legend('Location','northwest','box','off')

% ── Figure 3: TEMPORAL-STRUCTURE sensitivity — shuffle gap vs true τ ──────
if n_shuffle_recovery > 0
    gaps = [R.gap];
    figure('Color','w','Position',[80 80 560 470]); hold on
    plot(taus(finite), gaps(finite), 'o-','LineWidth',1.6)
    if any(~finite), yline(gaps(~finite), ':', 'accum', 'Color',[0.4 0.4 0.4]); end
    set(gca,'XScale','log'); grid on
    xlabel('True leak \tau (s)')
    ylabel('Shuffle gap:  CV-LL(real) − CV-LL(shuffled)')
    title('Sensitivity to temporal order vs integration horizon')
end

% ── Figure 4: distribution-match diagnostic (truncated CDFs vs real) ──────
figure('Color','w','Position',[80 80 700 470]); hold on
[fr,xr] = ecdf(real_durs); stairs(xr/fps, fr, 'k', 'LineWidth', 2.5);
for m = 1:numel(R)
    if numel(R(m).simdur) > 50
        [fs,xs] = ecdf(R(m).simdur); stairs(xs/fps, fs, 'Color', cmap(m,:), 'LineWidth', 1);
    end
end
xlabel('Freeze duration (s)'); ylabel('CDF'); xlim([0 maxT_fr/fps])
title('Truncated duration distributions: real (black) vs each τ (KS in console)')

% ── Figure 5: kernel-shape control — autocorrelation vs directed structure ─
% A recovered kernel that simply tracks the (symmetric) sm autocorrelation has no
% directed/causal content beyond "sm elevated around escape". ρ(kernel, sm-ACF)
% near 1 across τ ⇒ the hazard kernel is essentially the sm ACF, and the directed
% signal lives in the shuffle gap (Fig 3), not in the kernel shape.
rho_acf_all = arrayfun(@(k) R(k).res.rho_acf, 1:numel(R));
rho_eta_all = arrayfun(@(k) R(k).res.rho_eta, 1:numel(R));
figure('Color','w','Position',[80 80 980 430]);
tiledlayout(1,2,'TileSpacing','compact','Padding','compact')

nexttile; hold on
plot(taus(finite), rho_acf_all(finite), 'o-','LineWidth',1.6,'DisplayName','\rho(kernel, sm-ACF)')
plot(taus(finite), rho_eta_all(finite), 's--','LineWidth',1.2,'DisplayName','\rho(kernel, event excess)')
if any(~finite), yline(rho_acf_all(~finite), ':', 'accum', 'Color',[0.4 0.4 0.4]); end
set(gca,'XScale','log'); grid on; ylim([-0.2 1.05])
xlabel('True leak \tau (s)'); ylabel('shape correlation')
title('Recovered kernel ≈ sm autocorrelation?'); legend('Location','southwest','box','off')

nexttile; hold on   % shape overlay for the accumulator (or first τ): kernel vs ACF / event excess
kref = find(~finite,1); if isempty(kref), kref = 1; end
rr = R(kref).res; nrm = @(v) v / max(abs(v));
plot(rr.t_rel, nrm(rr.kernel),  'Color',[0.2 0.5 0.8],'LineWidth',2,  'DisplayName','recovered kernel')
plot(rr.t_rel, nrm(rr.acf_sm),  'k--','LineWidth',1.6,                 'DisplayName','sm autocorr R(\tau)')
plot(rr.t_rel, nrm(rr.eta_exc), 'Color',[0.85 0.4 0.2],'LineWidth',1.4,'DisplayName','event-triggered excess')
yline(0,'k:'); xline(0,'--k','escape')
xlabel('Time rel. to un-freeze (s)'); ylabel('shape (unit peak)')
if isinf(R(kref).tau), ttl='accumulator (τ=∞)'; else, ttl=sprintf('τ=%.2g s',R(kref).tau); end
title(sprintf('Shape overlay — %s', ttl)); legend('Location','best','box','off')

% ── Save ──────────────────────────────────────────────────────────────────
stamp = datestr(now,'yyyymmdd_HHMMSS'); %#ok<TNOW1,DATST>
try
    save(fullfile(paths.dataset, ['ground_truth_continuum_' stamp '.mat']), ...
        'R','tau_list','P','opts','real_med','real_durs','-v7.3');
    fprintf('\nSaved continuum recovery (stamp %s)\n', stamp);
catch ME
    warning('save failed: %s', ME.message);
end

% ════════════════════════════════════════════════════════════════════════
% Local functions
% ════════════════════════════════════════════════════════════════════════
function [dur, cens] = sim_leaky_all(theta, lambda, sm_in, P)
% One synthetic duration per bout under sim_leaky_accumulator(λ).
    nB = numel(sm_in); dur = nan(nB,1); cens = false(nB,1);
    for b = 1:nB
        drive = P.mu0 + P.beta * sm_in{b};
        k = sim_leaky_accumulator(drive, theta, P.sigma, lambda, P.dt, P.seed + b);
        if isnan(k), dur(b)=NaN; cens(b)=true; else, dur(b)=k; end
    end
end

function theta = calibrate_theta(simfun_of_theta, target_med, trunc_point)
% Bisection on θ so the median crossing duration (>trunc) ≈ target.
% Higher θ → harder to cross → longer durations / more censoring (monotone).
    lo = 0.02; hi = 80; theta = NaN;
    medfun = @(th) local_med(simfun_of_theta(th), trunc_point);
    for it = 1:30
        theta = 0.5*(lo+hi); m = medfun(theta);
        if isnan(m),            hi = theta;     % too few crossings → θ too high
        elseif m < target_med,  lo = theta;     % too fast → raise θ
        else,                   hi = theta;     % too slow → lower θ
        end
        if (hi-lo) < 1e-3, break; end
    end
end

function m = local_med(dur, trunc_point)
    valid = ~isnan(dur) & dur > trunc_point;
    if nnz(valid) < 50, m = NaN; else, m = median(dur(valid)); end
end
