clearvars; clc;
% ════════════════════════════════════════════════════════════════════════
% DDM KERNEL RECOVERY
%
% Simulate freeze durations from a leaky accumulator whose drift is a PARAMETRIC
% function of the real social motion,
%
%       drift(t) = mu0 + beta * sm(t),      leak lambda = 1/tau
%
% then recover the hazard kernel beta(tau) four ways and compare:
%
%   U  raised cosines, NO regularization        (reference / variance baseline)
%   A  raised cosines + L2 curvature penalty    (smoothness)
%   B  exponential bank + L1                    (sparse; selected atom = the leak)
%   C  raised cosines + L1                      (sparse, shape-neutral comparator)
%
% Sweeps the social gain BETA (include 0 for the null: sm does not drive the
% accumulator, so every kernel must collapse to ~0).
%
% Reuses, unchanged: generate_accumulator_dataset (DDM) -> build_hazard_design
% (person-period design) -> fit_kernels_sparse (the estimator used on real data,
% same code path as hazard_kernel_sparse.m).
% ════════════════════════════════════════════════════════════════════════

here     = fileparts(mfilename('fullpath'));
code_dir = fileparts(fileparts(fileparts(here)));                 % .../code
addpath(here);                                                    % generator, design, fitter
addpath(fullfile(code_dir, 'sims', 'simulators'));                % sim_leaky_accumulator
addpath(fullfile(code_dir, 'social_motion', 'hazard_kernel'));    % path_generator etc.

% ── kernel window / discretisation ──────────────────────────────────────────
fps             = 60;
kernel_past_s   = 2.0;    % causal history. Was 3.0. The measured true hazard kernel
                          % is at 1-2% of peak beyond 1 s and only ~10% of beta=4
                          % bouts have 3 s of in-freeze history at all, so the extra
                          % 2 s carried no signal and cost conditioning.
kernel_future_s = 0;    % acausal control (should recover ~0)
dt_frames       = 6;      % hazard sampled every 0.1 s
trunc_point     = 30;     % 0.5 s left truncation (generator AND design)
save_out        = true;
do_fit          = true;   % false -> simulate + draw the GENERATOR DIAGNOSTICS only
                          % (durations, sm-aligned-to-offset) and skip the ~1 h of
                          % fit_kernels_sparse. Those two figures depend on the
                          % generator alone, so nothing about them needs the fits.
rng(7);

% ── generative DDM: drift = mu0 + beta*sm, leak 1/tau ───────────────────────
P.dt          = 1/fps;
P.mu0         = 0.75;     % baseline urge to escape. MUST be > 0: at mu0 = 0 the
                          % beta = 0 null has exactly zero drift, so only diffusion
                          % ends a freeze and 66% of bouts never cross the bound --
                          % a degenerate null. Measured at theta = 3, n_aug = 1
                          % (cens/median_s at beta = 4/2/0, real median 2.15 s):
                          %   mu0=0.25  .02/.04/.38   1.68/2.58/4.77
                          %   mu0=0.50  .01/.01/.14   1.50/2.18/4.08
                          %   mu0=0.75  .00/.00/.03   1.35/1.87/3.30   <- chosen
                          %   mu0=1.00  .00/.00/.01   1.22/1.62/2.68
                          % 0.75 clears the censoring cliff with margin and the
                          % sweep still brackets the real median.
P.sigma       = 1.0;      % diffusion SD
P.theta       = 3;        % bound
P.seed        = 2305;
P.maxT_fr     = 1000;      % 10.5 s integration cap (no crossing -> censored)
P.trunc_point = trunc_point;
P.n_aug       = 2;        % DDM trajectories per real bout (more events)

beta_list = [4 2 0];   % PARAMETRIC social gain to sweep (0 = null)
tau_gen   = Inf;          % generative leak (s); Inf = perfect accumulator

% ── estimator options (mirror hazard_kernel_sparse.m) ───────────────────────
% use_nuisance = false -> design is [intercept, S*B] only: no baseline hazard
% b0(time-in-freeze), no mov/slo dummies. Duration dependence then has nowhere to
% go but the kernel, so expect a slow drift across lags on top of the true shape.
%
% The exponential bank is tau in [tau_lo, kernel_past_s] = [1/60, 1.0] s with 8
% atoms, replacing [0.05, 3.0] with 12. Measured on the reverse-correlation truth,
% cond(Bexp) drops 9.5e6 -> 3.0e3 (~3100x) and the best 2-atom fit improves
% 0.976 -> 0.998. Widening the bank UPWARD instead was strictly harmful: [0.05,10]
% degraded conditioning 55x for identical representational power, because over a
% short window the slow atoms are near-parallel straight lines. The gap was at the
% FAST end -- the true kernel is at 1/e by ~0.03 s, below the old shortest atom.
%
% nonneg_l1 constrains the L1 kernel coefficients (models B and C) to be >= 0, to
% stop near-duplicate atoms taking large opposite signs and cancelling. See the
% caveat in the fit_kernels_sparse header: under this constraint a non-negative
% kernel is NOT evidence that the kernel is positive.
opts = struct('nb_kernel',12, 'nb_acausal',0, 'bump_width',1.25, 'nb_baseline',6, ...
    'n_exp', 8, 'tau_lo', 1/fps, 'cv_folds', 5, 'kernel_past_s', kernel_past_s, ...
    'sel_rule', 'argmax', 'nonneg_l1', true, 'use_nuisance', false, 'verbose',true);

% ── load real bouts as TEMPLATE + motion cache ──────────────────────────────
% Only .fly/.onsets are used generatively (durations are re-simulated); the real
% sm supplies the drive, so regressor autocorrelation matches the real analysis.
threshold_imm = 2; threshold_mob = 2; threshold_pc = 4;
id_code = sprintf('imm%d_mob%d_pc%d', threshold_imm, threshold_mob, threshold_pc);
paths        = path_generator('folder','social_motion/hazard_kernel','bouts_id',id_code,'imfirst',false);
bouts        = importdata(fullfile(paths.dataset,  'bouts.mat'));
motion_cache = importdata(fullfile(paths.cache_path,'motion_cache.mat'));
thresholds   = define_thresholds('le_window', struct('le_window_sl',[5 55],'le_window_fl',[5 55]));
bouts        = bouts_formatting(bouts, thresholds);
bl_t = data_parser_new(bouts,'type','immobility','period','loom','window','le', ...
        'nloom',2:20,'min_dur',trunc_point);
nB = height(bl_t);
fprintf('Template: %d real freezes.\n', nB);

% ── global sm z-scoring (the scale the estimator reads back) ────────────────
flies = unique(bl_t.fly);  nn=0; s1=0; s2=0;
for i = 1:numel(flies)
    v = motion_cache(flies(i)); v = v(18000:end); v = v(~isnan(v));
    nn = nn + numel(v); s1 = s1 + sum(v); s2 = s2 + sum(v.^2);
end
sm_mu = s1/nn;  sm_sd = sqrt(s2/nn - sm_mu^2);

% ── per-bout DDM drive: real sm from onset forward, capped, raw (/10) units ──
% UNCENTRED: sm >= 0, so social motion only ever pushes TOWARD escape. The cost
% is that beta raises the mean drift as well as the kernel amplitude, so the
% panels of the sweep run at different speeds and are NOT rate-matched. Read the
% beta=0 panel as "no kernel AND slower process", not "no kernel" alone.
sm_in = cell(nB,1);
for b = 1:nB
    sm  = motion_cache(bl_t.fly(b)); sm = sm(:) ./ 10;
    on  = bl_t.onsets(b);
    cap = min(P.maxT_fr, numel(sm) - on + 1);
    sm_in{b} = sm(on : on + cap - 1);
    sm_in{b}(isnan(sm_in{b})) = 0;      % missing sm -> baseline drive only
end

% ── lag grid + approximate ground-truth shape ───────────────────────────────
back_fr = round(kernel_past_s*fps);  fwd_fr = round(kernel_future_s*fps);
lag_fr  = (-fwd_fr:back_fr)';  lag_s = lag_fr/fps;  caus = lag_s >= 0;
% A leaky accumulator's first-passage hazard weights history ~exp(-tau/tau_gen)
% (flat for tau_gen = Inf). This is the SHAPE reference only -- the bound and
% diffusion noise make the mapping approximate, so amplitude is fit by least
% squares before plotting and a small bias is expected.
gt = zeros(numel(lag_fr),1);
if isinf(tau_gen), gt(caus) = 1; else, gt(caus) = exp(-lag_s(caus)/tau_gen); end

cfg = struct('fps',fps, 'dt_frames',dt_frames, 'entry_fr',trunc_point, ...
    'lag_fr',lag_fr, 'sm_mu',sm_mu, 'sm_sd',sm_sd, 'grid_anchor','entry', 'mask_preonset', false);

% ── sweep beta: simulate -> design -> fit ───────────────────────────────────
R = struct('beta',{}, 'res',{}, 'info',{});
bl_all = cell(numel(beta_list),1);   % keep the synthetic bout tables for the
                                     % generator diagnostics below
fprintf('\n%-8s %8s %8s %8s %8s %8s\n', 'beta', 'events', 'q_U', 'q_A', 'q_B', 'q_C');
for m = 1:numel(beta_list)
    Pm = P;  Pm.beta = beta_list(m);
    Pm.lambda = 0;  if ~isinf(tau_gen), Pm.lambda = 1/tau_gen; end

    [bl, info] = generate_accumulator_dataset(bl_t, sm_in, Pm);
    bl_all{m}  = bl;
    if ~do_fit
        fprintf('%-8.3g %8d   (generation only — do_fit = false)\n', Pm.beta, height(bl));
        continue
    end

    [S, yev, tinf, mov, slo, gid] = build_hazard_design(bl, motion_cache, cfg);
    res = fit_kernels_sparse(S, yev, tinf, mov, slo, gid, lag_fr, fps, opts);

    rU = shape_fit(res.kernel_U, gt, caus);
    rA = shape_fit(res.kernel_A, gt, caus);
    rB = shape_fit(res.kernel_B, gt, caus);
    rC = shape_fit(res.kernel_C, gt, caus);
    fprintf('%-8.3g %8d %8.2f %8.2f %8.2f %8.2f\n', Pm.beta, res.n_events, rU, rA, rB, rC);
    if ~isempty(res.sel_tau)
        fprintf('         B selected tau atoms: %s s\n', num2str(res.sel_tau, '%.2f '));
    end

    R(m).beta = Pm.beta;  R(m).res = res;  R(m).info = info;
end

% ── figure: one panel per beta, all four estimators vs the reference shape ──
if do_fit
fh = figure('Color','w','Position',[60 60 450*numel(R) 400]);
tiledlayout(1, numel(R), 'TileSpacing','compact','Padding','compact');
for m = 1:numel(R)
    rr = R(m).res;  nexttile; hold on
    apply_generic(gca, 'ylim', [-0.05 0.25], 'xlim', [-2.05 0.05])
    % shade the acausal control region (should sit at ~0)
    % Keep explicit handles: the patch below is uistack'd to the bottom of the
    % child list, and a legend given only labels would hand it the first one and
    % shift every curve's label by one.
    h(1) = plot(rr.t_rel, rr.kernel_U, '-',  'Color',[0.6 0.6 0.6], 'LineWidth',1);
    h(2) = plot(rr.t_rel, rr.kernel_A, '-',  'Color',[0.2 0.5 0.8], 'LineWidth',2.3);
    h(3) = plot(rr.t_rel, rr.kernel_B, '-',  'Color',[0.85 0.4 0.2], 'LineWidth',2.3);
    h(4) = plot(rr.t_rel, rr.kernel_C, '-', 'Color',[0.3 0.3 0.3], 'LineWidth',2.3);
    yl = ylim; tac = max(rr.t_rel);
    hp = patch([0 tac tac 0],[yl(1) yl(1) yl(2) yl(2)],[.93 .93 .93], ...
        'EdgeColor','none', 'HandleVisibility','off');
    uistack(hp,'bottom'); ylim(yl);
    yline(0,'k--', 'LineWidth', 1); xline(0,'--k', 'LineWidth', 1);
    xlabel('Time to unfreeze (s)');
    if m==1
        ylabel('\beta(\tau)');
        legend(h, {'A','A-L2','B-L1','C-L1'}, ...
            'Location','northwest','box','off', 'FontSize', 18);
    end
    title(sprintf('\\beta_{gen}=%.3g', R(m).beta), 'FontSize', 20);
end

% ── figure: the weighted basis components b_j*B(:,j) that sum to A and B ────
% Works standalone off any saved run, so graphics can be tweaked without refitting:
%   plot_kernel_components('results/ddm_kernel_recovery_<stamp>.mat')
fh_cmp = plot_kernel_components(R, opts, lag_fr, fps, beta_list);
end

% ── generator diagnostics: what was simulated, and what the design sees ─────
% Both depend on the GENERATOR only (beta, mu0, theta), never on mask_preonset /
% nb_baseline / use_nuisance -- so they are identical across estimator variants
% and the panels vary with beta alone.
[fh_dur, fh_sm] = plot_generator_diagnostics(bl_all, bl_t, motion_cache, beta_list, cfg, P);

if save_out
    stamp  = datestr(now,'yyyymmdd_HHMMSS'); %#ok<TNOW1,DATST>
    outdir = fullfile(here,'results');
    if ~exist(outdir,'dir'), mkdir(outdir); end
    exportgraphics(fh_dur, fullfile(outdir, ['ddm_recovery_durations_' stamp '.png']), 'Resolution',200);
    exportgraphics(fh_sm,  fullfile(outdir, ['ddm_recovery_smalign_'   stamp '.png']), 'Resolution',200);
    if do_fit    % no kernels to save when only the generator diagnostics were drawn
        save(fullfile(outdir, ['ddm_kernel_recovery_' stamp '.mat']), ...
            'R','P','beta_list','tau_gen','opts','cfg','lag_fr','gt','bl_all','-v7.3');
        exportgraphics(fh, fullfile(outdir, ['ddm_kernel_recovery_' stamp '.png']), 'Resolution',200);
        exportgraphics(fh_cmp, fullfile(outdir, ['ddm_recovery_components_' stamp '.png']), 'Resolution',200);
    end
    fprintf('\nSaved to %s (stamp %s)\n', outdir, stamp);
end

% ════════════════════════════════════════════════════════════════════════
function q = shape_fit(k, gt, caus)
% Agreement with the ground-truth SHAPE over causal lags, amplitude-free:
%   q = 1 - ||k - a*gt|| / ||k||,   a = least-squares scale (as in scale_gt)
% 1 = perfect, 0 = no better than nothing. A plain correlation cannot be used
% here: for tau_gen = Inf the truth is CONSTANT over causal lags, so corr(k,gt)
% is 0/0 = NaN for every model. This reduces to "how flat is k" in that case and
% to a shape match for finite tau_gen, so one column reads the same either way.
    kc = k(caus);  gc = gt(caus);
    if norm(kc) < eps, q = NaN; return; end
    a = (gc' * kc) / max(gc' * gc, eps);
    q = 1 - norm(kc - a*gc) / norm(kc);
end

function g = scale_gt(gt, k, caus)
% Least-squares scale the ground-truth shape onto a recovered kernel, so the
% comparison is about SHAPE (the DDM->kernel amplitude map is not identity).
    a = (gt(caus)' * k(caus)) / max(gt(caus)' * gt(caus), eps);
    g = a * gt;
end
