clearvars; close all; clc;
% ════════════════════════════════════════════════════════════════════════
% KERNEL RECOVERY — sparse / exponential methods
%
% Control experiment for hazard_kernel_sparse.m. We SIMULATE un-freeze events
% under a KNOWN hazard kernel β(τ) (simulate_hazard_bouts), driven by the REAL
% social motion each fly saw, then recover the kernel with the SAME estimator the
% real analysis uses (fit_kernels_sparse). This is the paper's own validation
% strategy (Mineault, Barthelmé & Pack 2009): put a template in, check it comes
% back out, and measure how FEW events are needed.
%
% Two experiments:
%   1. Shape recovery — integrator / transient / null ground truths. Does the
%      sparse exponential model (B) recover the true shape (corr with truth),
%      the true leak τ, and correctly return ≈0 for the null (specificity)?
%   2. Efficiency — subsample bouts to vary event count and plot recovery corr
%      vs #events for L2 smoothness (A) vs L1 exponential (B) vs L1 flexible (C).
%      The paper's claim is that the sparse model recovers with fewer events.
% ════════════════════════════════════════════════════════════════════════

here = fileparts(mfilename('fullpath'));
addpath(here);   % fit_kernels_sparse, simulate_hazard_bouts

fps = 60;
kernel_past_s   = 3.0;
kernel_future_s = 1.0;
dt_frames   = 6;
trunc_point = 30;   entry_fr = trunc_point;
save_out    = true;
rng(7);

% ── load real bouts + motion cache (same pipeline as hazard_kernel_sparse.m) ──
threshold_imm = 2; threshold_mob = 2; threshold_pc = 4;
id_code = sprintf('imm%d_mob%d_pc%d', threshold_imm, threshold_mob, threshold_pc);
paths        = path_generator('folder','social_motion/hazard_kernel','bouts_id',id_code,'imfirst',false);
bouts        = importdata(fullfile(paths.dataset,  'bouts.mat'));
motion_cache = importdata(fullfile(paths.cache_path,'motion_cache.mat'));
thresholds = define_thresholds('le_window', struct('le_window_sl',[5 55],'le_window_fl',[5 55]));
bouts = bouts_formatting(bouts, thresholds);
bl_t  = data_parser_new(bouts, 'type','immobility','period','loom','window','le','nloom',2:20,'min_dur',trunc_point);
bl_t  = impose_contact_threshold(bl_t, 'threshold', 75, 'type','onlyfreeze');
bl_t.dur_frames = bl_t.ending_time;            % follow-up cap per bout
bl_t = bl_t(bl_t.dur_frames > trunc_point, :);
nB = height(bl_t);
fprintf('Template: %d real freezes.\n', nB);

% ── global z-scoring of sm (as in the real analysis) ─────────────────────────
flies = unique(bl_t.fly);  nn=0; s1=0; s2=0;
for i = 1:numel(flies)
    v = motion_cache(flies(i)); v = v(:); v = v(~isnan(v));
    nn = nn + numel(v); s1 = s1 + sum(v); s2 = s2 + sum(v.^2);
end
sm_mu = s1/nn;  sm_sd = sqrt(s2/nn - sm_mu^2);

% ── lag grid + ground-truth kernels (causal only; acausal truth = 0) ─────────
back_fr = round(kernel_past_s*fps);  fwd_fr = round(kernel_future_s*fps);
lag_fr  = (-fwd_fr:back_fr)';  lag_s = lag_fr/fps;  caus = lag_s >= 0;  nLag = numel(lag_fr);

truths = struct('name',{}, 'beta',{}, 'tau',{});
bi = zeros(nLag,1);  bi(caus) = exp(-lag_s(caus)/0.80);            % integrator, τ=0.8 s
truths(1) = struct('name','integrator (τ=0.80s)', 'beta',bi, 'tau',0.80);
bt = zeros(nLag,1);  bt(caus) = exp(-lag_s(caus)/0.06);            % transient, τ=0.06 s
truths(2) = struct('name','transient (τ=0.06s)',  'beta',bt, 'tau',0.06);
bn = zeros(nLag,1);                                                % null
truths(3) = struct('name','null (β=0)',           'beta',bn, 'tau',NaN);

% ── simulation + estimation options ──────────────────────────────────────────
sp = struct('dt_frames',dt_frames, 'entry_fr',entry_fr, 'n_aug',1, 'seed',100, ...
    'target_eta_sd',0.6, 'target_hazard',0.02, 'sm_mu',sm_mu, 'sm_sd',sm_sd, ...
    'mask_preonset',false, ...   % true → zero out lags reaching before freeze onset
    'grid_anchor','entry');      % 'entry' from onset | 'offset' from freeze break
opts = struct('nb_kernel',12,'nb_acausal',4,'bump_width',1.25,'nb_baseline',6, ...
    'n_exp',10,'cv_folds',5,'kernel_past_s',kernel_past_s,'sel_rule','argmax','verbose',false);

% ════════════════════════════════════════════════════════════════════════
% Experiment 1 — shape recovery for each ground-truth kernel
% ════════════════════════════════════════════════════════════════════════
fprintf('\n══ Experiment 1: shape recovery ══\n');
fprintf('%-22s %8s %8s %8s %8s\n', 'truth', 'events', 'r_A', 'r_B(exp)', 'r_C');
fh1 = figure('Color','w','Position',[80 80 1200 380]);
tiledlayout(1, numel(truths), 'TileSpacing','compact','Padding','compact');
for ti = 1:numel(truths)
    spi = sp;  spi.seed = 100 + ti;
    [S,yev,tinf,mov,slo,gid,info] = simulate_hazard_bouts(bl_t, motion_cache, truths(ti).beta, lag_fr, fps, spi);
    R = fit_kernels_sparse(S, yev, tinf, mov, slo, gid, lag_fr, fps, opts);

    bt_c = info.beta_true(caus);
    [rA, tauA] = recov_metrics(R.kernel_A, bt_c, caus, lag_s);
    [rB, tauB] = recov_metrics(R.kernel_B, bt_c, caus, lag_s);
    [rC, tauC] = recov_metrics(R.kernel_C, bt_c, caus, lag_s);

    fprintf('%-22s %8d %8.2f %8.2f %8.2f', truths(ti).name, info.n_events, rA, rB, rC);
    if ~isnan(truths(ti).tau)
        fprintf('   | true τ=%.2f  recovered τ_com: A=%.2f B=%.2f C=%.2f', truths(ti).tau, tauA, tauB, tauC);
    else
        fprintf('   | null: max|β̂|  A=%.3f B=%.3f C=%.3f', max(abs(R.kernel_A(caus))), ...
            max(abs(R.kernel_B(caus))), max(abs(R.kernel_C(caus))));
    end
    if ~isempty(R.sel_tau), fprintf('  [B atoms: %s]', num2str(R.sel_tau, '%.2f ')); end
    fprintf('\n');

    nexttile; hold on
    plot(R.t_rel, info.beta_true, 'k-', 'LineWidth', 2.5);
    plot(R.t_rel, R.kernel_A, 'Color',[0.2 0.5 0.8], 'LineWidth', 1.6);
    plot(R.t_rel, R.kernel_B, 'Color',[0.85 0.4 0.2], 'LineWidth', 1.8);
    plot(R.t_rel, R.kernel_C, '--', 'Color',[0.3 0.3 0.3], 'LineWidth', 1.2);
    yline(0,'k:'); xline(0,'--k');
    xlabel('time rel. to un-freeze (s)'); ylabel('\beta(\tau)');
    title(sprintf('%s  (n=%d ev)', truths(ti).name, info.n_events));
    if ti==1, legend({'truth','A: L2','B: L1 exp','C: L1 flex'}, 'Location','best','box','off'); end
end

% ════════════════════════════════════════════════════════════════════════
% Experiment 2 — recovery vs #events (subsample bouts): the efficiency curve
% ════════════════════════════════════════════════════════════════════════
fprintf('\n══ Experiment 2: recovery vs #events (integrator truth) ══\n');
frac_list = [0.1 0.2 0.4 0.7 1.0];
nfr = numel(frac_list);
ev_n = nan(nfr,1);  rr_A = nan(nfr,1);  rr_B = nan(nfr,1);  rr_C = nan(nfr,1);
for fi = 1:nfr
    idx = randperm(nB, max(20, round(frac_list(fi)*nB)));
    bl_s = bl_t(idx,:);
    spi = sp;  spi.seed = 500 + fi;
    [S,yev,tinf,mov,slo,gid,info] = simulate_hazard_bouts(bl_s, motion_cache, truths(1).beta, lag_fr, fps, spi);
    R = fit_kernels_sparse(S, yev, tinf, mov, slo, gid, lag_fr, fps, opts);
    bt_c = info.beta_true(caus);
    ev_n(fi) = info.n_events;
    rr_A(fi) = recov_metrics(R.kernel_A, bt_c, caus, lag_s);
    rr_B(fi) = recov_metrics(R.kernel_B, bt_c, caus, lag_s);
    rr_C(fi) = recov_metrics(R.kernel_C, bt_c, caus, lag_s);
    fprintf('frac=%.2f  events=%5d  r: A=%.2f  B(exp)=%.2f  C=%.2f\n', frac_list(fi), ev_n(fi), rr_A(fi), rr_B(fi), rr_C(fi));
end

[ev_s, ord] = sort(ev_n);
fh2 = figure('Color','w','Position',[100 100 640 460]); hold on
plot(ev_s, rr_A(ord), 'o-', 'Color',[0.2 0.5 0.8], 'LineWidth',2, 'MarkerFaceColor','w');
plot(ev_s, rr_B(ord), 's-', 'Color',[0.85 0.4 0.2], 'LineWidth',2, 'MarkerFaceColor','w');
plot(ev_s, rr_C(ord), '^--','Color',[0.3 0.3 0.3], 'LineWidth',1.6, 'MarkerFaceColor','w');
ylim([0 1]); grid on
xlabel('number of un-freeze events'); ylabel('recovery corr  r(\betâ, \beta_{true})');
title('Kernel recovery vs #events (integrator truth)');
legend({'A: L2 smoothness','B: L1 exponential','C: L1 flexible'}, 'Location','southeast','box','off');

% ── persist ───────────────────────────────────────────────────────────────
if save_out
    stamp  = datestr(now, 'yyyymmdd_HHMMSS'); %#ok<TNOW1,DATST>
    outdir = fullfile(here, 'results');
    if ~exist(outdir,'dir'), mkdir(outdir); end
    try
        save(fullfile(outdir, ['recovery_sparse_' stamp '.mat']), ...
            'truths','frac_list','ev_n','rr_A','rr_B','rr_C','sp','opts', ...
            'kernel_past_s','kernel_future_s','dt_frames', '-v7.3');
        exportgraphics(fh1, fullfile(outdir, ['recovery_shapes_'    stamp '.png']), 'Resolution',180);
        exportgraphics(fh2, fullfile(outdir, ['recovery_efficiency_' stamp '.png']), 'Resolution',180);
        fprintf('\nSaved recovery results + figures to %s (stamp %s)\n', outdir, stamp);
    catch ME
        warning('save_out failed: %s', ME.message);
    end
end

% ════════════════════════════════════════════════════════════════════════
function [r, tau_com] = recov_metrics(kernel, beta_true_caus, caus, lag_s)
% Pearson corr between recovered and true CAUSAL kernel (NaN if truth is flat),
% and the recovered leak as the center-of-mass lag of the positive causal weight.
    k_c = kernel(caus);
    if std(beta_true_caus) < eps
        r = NaN;
    else
        a = k_c - mean(k_c);  b = beta_true_caus - mean(beta_true_caus);
        r = (a'*b) / (norm(a)*norm(b) + eps);
    end
    kpos = max(k_c, 0);  ls = lag_s(caus);
    if sum(kpos) > 0, tau_com = sum(ls .* kpos) / sum(kpos); else, tau_com = NaN; end
end
