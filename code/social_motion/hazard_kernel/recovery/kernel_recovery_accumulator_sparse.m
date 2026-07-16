clearvars; close all; clc;
% ════════════════════════════════════════════════════════════════════════
% BOUNDED-ACCUMULATOR RECOVERY — sparse / exponential estimator
%
% The mechanistic control. Freeze durations are the FIRST-PASSAGE times of a
% genuine bounded leaky accumulator (sim_leaky_accumulator) integrating the REAL
% social motion each fly saw:
%       x(k+1) = x(k)(1 - λΔt) + (μ0 + β·sm(k))Δt + σ√Δt·randn,   λ = 1/τ
%       un-freeze = first k with x ≥ θ.
% We sweep the generative leak time-constant τ, run the SAME sparse-exponential
% estimator used on real data (build_hazard_design → fit_kernels_sparse), and ask
% whether the recovered kernel is exponential with the RIGHT τ. Unlike the
% kernel-level control (simulate_hazard_bouts), the bound + diffusion noise +
% first-passage make the τ→kernel mapping only approximate — so recovering τ here
% is a real test of the biological interpretation, not a perfectly-specified fit.
%
% Reuses your existing generator (generate_accumulator_dataset + sim_leaky_
% accumulator); only the ESTIMATOR is swapped to the new sparse/exponential one.
% ════════════════════════════════════════════════════════════════════════

here     = fileparts(mfilename('fullpath'));
code_dir = fileparts(fileparts(fileparts(here)));           % .../code
addpath(here);                                              % fit_kernels_sparse, build_hazard_design, generate_accumulator_dataset
addpath(fullfile(code_dir, 'sims', 'simulators'));          % sim_leaky_accumulator

fps = 60;
kernel_past_s   = 3.0;
kernel_future_s = 1.0;
dt_frames   = 6;
trunc_point = 30;   entry_fr = trunc_point;
save_out    = true;
rng(11);

% ── generative accumulator regime (drive = μ0 + β·sm) ────────────────────────
P.dt          = 1/fps;
P.mu0         = 0;      % baseline drive -> duration timescale
P.beta        = 5;      % social gain
P.sigma       = 1.0;    % diffusion SD
P.theta       = 2;      % bound
P.n_aug       = 1;      % synthetic trajectories per real bout
P.maxT_fr     = 630;    % 10.5 s integration cap (no crossing -> censored)
P.trunc_point = trunc_point;
P.seed        = 24200;

% swept leak time-constants (s). Inf = perfect accumulator (flat kernel).
tau_list = [Inf 2.0 1.0 0.5 0.25];

% ── load real bouts + motion cache (same pipeline as hazard_kernel_sparse.m) ──
threshold_imm = 2; threshold_mob = 2; threshold_pc = 4;
id_code = sprintf('imm%d_mob%d_pc%d', threshold_imm, threshold_mob, threshold_pc);
paths        = path_generator('folder','social_motion/hazard_kernel','bouts_id',id_code,'imfirst',false);
bouts        = importdata(fullfile(paths.dataset,  'bouts.mat'));
motion_cache = importdata(fullfile(paths.cache_path,'motion_cache.mat'));
thresholds = define_thresholds('le_window', struct('le_window_sl',[5 55],'le_window_fl',[5 55]));
bouts = bouts_formatting(bouts, thresholds);
bl_t  = data_parser_new(bouts, 'type','immobility','period','loom','window','le','nloom',2:20,'min_dur',trunc_point);
nB = height(bl_t);
fprintf('Template: %d real freezes.\n', nB);

% ── per-bout sm drive (onset forward, capped, /10, NaN->0) — as in the sweep ──
sm_in = cell(nB,1);
for b = 1:nB
    sm  = motion_cache(bl_t.fly(b)); sm = sm(:) ./ 10;
    on  = bl_t.onsets(b);
    cap = min(P.maxT_fr, numel(sm) - on + 1);
    sm_in{b} = sm(on : on + cap - 1);
    sm_in{b}(isnan(sm_in{b})) = 0;
end

% ── global z-scoring of raw sm (for the ESTIMATOR's design) ──────────────────
flies = unique(bl_t.fly);  nn=0; s1=0; s2=0;
for i = 1:numel(flies)
    v = motion_cache(flies(i)); v = v(:); v = v(~isnan(v));
    nn = nn + numel(v); s1 = s1 + sum(v); s2 = s2 + sum(v.^2);
end
sm_mu = s1/nn;  sm_sd = sqrt(s2/nn - sm_mu^2);

% ── lag grid + estimator config ──────────────────────────────────────────────
back_fr = round(kernel_past_s*fps);  fwd_fr = round(kernel_future_s*fps);
lag_fr  = (-fwd_fr:back_fr)';  lag_s = lag_fr/fps;  caus = lag_s >= 0;
cfg  = struct('fps',fps,'dt_frames',dt_frames,'entry_fr',entry_fr,'lag_fr',lag_fr,'sm_mu',sm_mu,'sm_sd',sm_sd, ...
    'mask_preonset',false, ...   % true → zero out lags reaching before freeze onset
    'grid_anchor','entry');      % 'entry' from onset | 'offset' from freeze break
opts = struct('nb_kernel',12,'nb_acausal',4,'bump_width',1.25,'nb_baseline',6, ...
    'n_exp',10,'cv_folds',5,'kernel_past_s',kernel_past_s,'sel_rule','argmax','verbose',false);

% ════════════════════════════════════════════════════════════════════════
% Sweep τ: generate bounded-accumulator freezes -> recover with sparse estimator
% ════════════════════════════════════════════════════════════════════════
nT = numel(tau_list);
gen_tau = nan(nT,1); rec_tauhat = nan(nT,1); rec_taucom = nan(nT,1);
R2 = nan(nT,1); nev = nan(nT,1); kernels = nan(numel(lag_fr), nT);
fprintf('\n%-14s %8s %10s %10s %8s\n','gen τ','events','τ̂_exp(B)','τ_com(B)','R²');
for m = 1:nT
    tau = tau_list(m);
    Pm  = P;  Pm.lambda = 0;  if ~isinf(tau), Pm.lambda = 1/tau; end
    Pm.seed = P.seed + m*777;

    bl = generate_accumulator_dataset(bl_t, sm_in, Pm);
    [S,yev,tinf,mov,slo,gid] = build_hazard_design(bl, motion_cache, cfg);
    Rk = fit_kernels_sparse(S, yev, tinf, mov, slo, gid, lag_fr, fps, opts);

    kernels(:,m) = Rk.kernel_B;
    [tauh, r2]   = fit_leak_tau(Rk.kernel_B, lag_s, caus, Rk.kernel_B_se);
    taucom       = com_tau(Rk.kernel_B, lag_s, caus);
    gen_tau(m)=tau; rec_tauhat(m)=tauh; rec_taucom(m)=taucom; R2(m)=r2; nev(m)=Rk.n_events;

    fprintf('%-14s %8d %10.2f %10.2f %8.2f', ternary(isinf(tau),'Inf',sprintf('%.2f',tau)), ...
        Rk.n_events, tauh, taucom, r2);
    if ~isempty(Rk.sel_tau), fprintf('   [B atoms: %s]', num2str(Rk.sel_tau,'%.2f ')); end
    fprintf('\n');
end

% ── β = 0 baseline (social motion causally irrelevant): the null floor ───────
% Drive = μ0 only, so un-freeze times are INDEPENDENT of sm. If the recovered
% kernel here is ≈0, the τ-sweep kernels above are genuine social-motion effects;
% if it still shows a ~0.45 s kernel, that timescale is a bound/first-passage
% artifact (sm autocorrelation leaking in), not a causal social-motion effect.
Pm0 = P;  Pm0.beta = 0;  Pm0.lambda = 0;  Pm0.seed = P.seed + 99;
bl0 = generate_accumulator_dataset(bl_t, sm_in, Pm0);
[S0,y0,ti0,mv0,sl0,g0] = build_hazard_design(bl0, motion_cache, cfg);
R0 = fit_kernels_sparse(S0, y0, ti0, mv0, sl0, g0, lag_fr, fps, opts);
signal_peak = max(abs(kernels(:)), [], 'omitnan');
fprintf('\n── β=0 null (social motion causally irrelevant) ──\n');
fprintf('events=%d   max|β̂| causal:  A=%.3f  B=%.3f  C=%.3f   (signal peak across sweep ~%.2f)\n', ...
    R0.n_events, max(abs(R0.kernel_A(caus))), max(abs(R0.kernel_B(caus))), ...
    max(abs(R0.kernel_C(caus))), signal_peak);
if ~isempty(R0.sel_tau)
    fprintf('   B atoms: %s   (a lone boundary atom with ~0 weight is the benign catch-all)\n', num2str(R0.sel_tau,'%.2f '));
else
    fprintf('   B selected no atoms (kernel ≈ 0) — clean null.\n');
end
kernel0_B = R0.kernel_B;   %#ok<NASGU>  saved below

% ── Fig 1: recovered kernels (should broaden as τ grows) ─────────────────────
fh1 = figure('Color','w','Position',[80 80 720 460]); hold on
cmap = parula(nT);
lg = cell(nT,1);
for m = 1:nT
    plot(-lag_s, kernels(:,m), 'Color', cmap(m,:), 'LineWidth', 2);
    lg{m} = ternary(isinf(tau_list(m)), 'τ=Inf (flat)', sprintf('τ=%.2f s', tau_list(m)));
end
plot(-lag_s, R0.kernel_B, 'k:', 'LineWidth', 2);  lg{end+1} = '\beta=0 null';
yline(0,'k:'); xline(0,'--k','escape')
xlabel('time rel. to un-freeze (s)   —   past < 0'); ylabel('recovered \beta(\tau)  (model B)')
title('Recovered kernels vs generative leak τ (bounded accumulator)')
legend(lg, 'Location','best', 'box','off')

% ── Fig 3: the OTHER component — baseline hazard b0(time-in-freeze) ──────────
% From the β=0 null (no sm effect), so b0 carries all the timing: this is the
% "probability of breaking at time t" that the estimator separates from β(τ).
[ts, ord] = sort(R0.tinf);
fh3 = figure('Color','w','Position',[120 120 640 440]); hold on
plot(ts, 1 ./ (1 + exp(-R0.base_logit(ord))), 'Color',[0.2 0.5 0.8], 'LineWidth', 2);
xlabel('time in freeze (s)'); ylabel('baseline hazard  P(unfreeze \mid frozen, sm=0)')
title('Baseline hazard b_0(t): probability of breaking vs time in freeze (\beta=0 null)')
grid on

% ── Fig 2: recovered τ vs generative τ against the identity line ─────────────
fin = isfinite(gen_tau);
fh2 = figure('Color','w','Position',[100 100 560 480]); hold on
mx = max([gen_tau(fin); rec_tauhat(fin); 2.5]);
plot([0 mx],[0 mx],'k--','LineWidth',1);
plot(gen_tau(fin), rec_tauhat(fin), 'o-', 'Color',[0.85 0.4 0.2], 'LineWidth',2, 'MarkerFaceColor','w');
plot(gen_tau(fin), rec_taucom(fin), 's--','Color',[0.2 0.5 0.8], 'LineWidth',1.5, 'MarkerFaceColor','w');
axis equal; xlim([0 mx]); ylim([0 mx]); grid on
xlabel('generative leak \tau (s)'); ylabel('recovered \tau (s)')
title('Leak-time-constant recovery (bounded accumulator -> sparse estimator)')
legend({'identity','\taû exponential fit (B)','\tau_{com} (B)'}, 'Location','northwest','box','off')

% ── persist ───────────────────────────────────────────────────────────────
if save_out
    stamp  = datestr(now, 'yyyymmdd_HHMMSS'); %#ok<TNOW1,DATST>
    outdir = fullfile(here, 'results');
    if ~exist(outdir,'dir'), mkdir(outdir); end
    try
        save(fullfile(outdir, ['recovery_accum_sparse_' stamp '.mat']), ...
            'tau_list','gen_tau','rec_tauhat','rec_taucom','R2','nev','kernels', ...
            'lag_s','P','opts','kernel_past_s','kernel0_B', '-v7.3');
        exportgraphics(fh1, fullfile(outdir, ['recovery_accum_kernels_'  stamp '.png']), 'Resolution',180);
        exportgraphics(fh2, fullfile(outdir, ['recovery_accum_tau_'      stamp '.png']), 'Resolution',180);
        exportgraphics(fh3, fullfile(outdir, ['recovery_accum_baseline_' stamp '.png']), 'Resolution',180);
        fprintf('\nSaved accumulator-recovery results + figures to %s (stamp %s)\n', outdir, stamp);
    catch ME
        warning('save_out failed: %s', ME.message);
    end
end

% ════════════════════════════════════════════════════════════════════════
function [tau_hat, R2] = fit_leak_tau(kernel, lag_s, caus, kernel_se)
% Weighted exponential fit β(τ) ≈ A·exp(-τ/τ_leak)+c to the causal kernel; returns
% τ_leak (unbiased by the finite window, unlike center-of-mass). Mirrors the leak
% fit in hazard_kernel_sm.m.
    tau_s = lag_s(caus);  kc = kernel(caus);
    if nargin >= 4 && ~isempty(kernel_se) && any(kernel_se(caus) > 0)
        sec = kernel_se(caus);  w = 1 ./ max(sec, eps).^2;
    else
        w = ones(size(kc));
    end
    expmod = @(p,t) p(1).*exp(-t./exp(p(2))) + p(3);
    sse    = @(p) sum(w .* (kc - expmod(p,tau_s)).^2);
    p0     = [max(kc)-min(kc), log(0.5), min(kc)];
    opt    = optimset('Display','off','MaxFunEvals',1e4,'MaxIter',1e4);
    pf     = fminsearch(sse, p0, opt);
    tau_hat = exp(pf(2));
    res  = kc - expmod(pf, tau_s);  wmean = sum(w.*kc)/sum(w);
    R2   = 1 - sum(w.*res.^2) / sum(w.*(kc - wmean).^2);
end

function tc = com_tau(kernel, lag_s, caus)
    kc = max(kernel(caus), 0);  ls = lag_s(caus);
    if sum(kc) > 0, tc = sum(ls.*kc)/sum(kc); else, tc = NaN; end
end

function out = ternary(c, a, b)
    if c, out = a; else, out = b; end
end
