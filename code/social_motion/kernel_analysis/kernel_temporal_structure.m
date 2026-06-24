clearvars; clc;
% ════════════════════════════════════════════════════════════════════════
% DOES THE TEMPORAL STRUCTURE OF SOCIAL MOTION MATTER?
%
% Hold the generative parameters FIXED (gain, leak, bound) and vary only how the
% drive is CONSTRUCTED, then run the identical hazard-kernel analysis
% (run_hazard_kernel) and compare the recovered kernels. Drive conditions:
%
%   'timeseries' — accumulator integrates the full real sm(t). Temporal
%                  structure is present and causal. Reference.
%   'mean'       — accumulator integrates a CONSTANT per-bout drift
%                  mu0 + beta*mean(sm_bout): the time series is collapsed to one
%                  scalar, so NO within-bout timing carries information. The
%                  analysis still reads the real time-resolved sm.
%   'white'      — drive = analysed signal = white noise (no autocorrelation).
%   'ar1'        — drive = analysed signal = AR(1) with autocorrelation time tau_ac.
%
% WHAT TO EXPECT
%   - 'timeseries' (esp. white): kernel is CAUSAL — weight in the past (tau<=0),
%     ~0 in the acausal/future region (tau>0).
%   - 'mean': kernel goes FLAT and SYMMETRIC across tau=0 — the acausal half is
%     as large as the causal half, because any sm sample (past OR future) is just
%     a proxy for the bout mean. Strong acausal weight = "only the level matters,
%     timing is irrelevant". That is the discriminator (NOT the time-shuffle: a
%     perfect integrator's sum is order-invariant, so its shuffle gap is ~0 too).
%
% Headline metric: acausal/causal mean |beta(tau)| ratio. ~0 => causal/time-
% resolved; ~1 => level-only / no temporal structure.
% ════════════════════════════════════════════════════════════════════════

% ── paths ──────────────────────────────────────────────────────────────────
here     = fileparts(mfilename('fullpath'));
code_dir = fileparts(fileparts(here));
addpath(fullfile(code_dir, 'sims', 'simulators'));            % sim_leaky_accumulator
addpath(fullfile(code_dir, 'social_motion', 'proximity'));    % run_hazard_kernel

% ── FIXED generative regime (only the drive CONSTRUCTION varies below) ──────
fps           = 60;
P.dt          = 1/fps;
P.mu0         = 0;
P.beta        = 8;       % social gain (fixed across conditions)
P.sigma       = 1.0;
P.seed        = 222;
P.maxT_fr     = 630;
P.trunc_point = 0;
P.n_aug       = 1;
P.theta       = 8.5;
tau_gen       = Inf;     % generative leak (Inf = perfect accumulator); fixed
P.lambda      = 0; if ~isinf(tau_gen), P.lambda = 1/tau_gen; end

% ── drive conditions to compare ─────────────────────────────────────────────
conditions = {'timeseries', 'mean'};   % add 'white','ar1' to also test autocorrelation
tau_ac     = 0.5;                       % AR(1) autocorrelation time (s) for 'ar1'

calibrate  = false;     % calibrate theta per condition to the real median
save_out   = true;
do_eta_overlay = true;  % extra figure: raw kernel vs event-triggered sm excess (eta_exc)

% ── hazard-analysis options (kernel_future_s>0 so the acausal control exists) ─
opts = struct('fps',fps, 'kernel_past_s', 3, 'kernel_future_s', 2, 'dt_frames', 3, ...
    'nb_kernel', 9, 'nb_acausal', 6, 'bump_width', 1.25, 'nb_baseline', 6, ...
    'cv_folds', 5, 'n_shuffle', 0, 'shuffle_mode', 'full', 'trunc_point', P.trunc_point, ...
    'do_plot', true, 'do_stage2', false, 'do_unpen', true, 'use_penalty', true, 'anchor_zero', false);

% ── load real bouts + motion cache (same pipeline as hazard_kernel_sm.m) ────
threshold_imm = 2; threshold_mob = 2; threshold_pc = 4;
id_code = sprintf('imm%d_mob%d_pc%d', threshold_imm, threshold_mob, threshold_pc);
paths = path_generator('folder','social_motion/proximity','bouts_id',id_code,'imfirst',false);
bouts = importdata(fullfile(paths.dataset,  'bouts.mat'));
motion_cache = importdata(fullfile(paths.cache_path, 'motion_cache.mat'));
thresholds = define_thresholds('le_window', struct('le_window_sl',[5 55],'le_window_fl',[5 55]));
bouts = bouts_formatting(bouts, thresholds);
bl_t  = data_parser_new(bouts, 'type','immobility', 'period','loom', 'window','le', ...
        'nloom', 2:20, 'min_dur', P.trunc_point);

if ismember('durations', bl_t.Properties.VariableNames)
    real_dur = double(bl_t.durations);
else
    real_dur = double(bl_t.ends - bl_t.onsets);
end

real_med  = median(real_dur);
fprintf('Template: %d real freezes, median %.0f fr (%.2f s)\n', height(bl_t), real_med, real_med/fps);

% ── reference scale: SD of the real (/10) sm over the analysed windows ──────
% Used to give the white/AR(1) drives a comparable amplitude to the real sm.
sd_ref = ref_sm_sd(motion_cache, bl_t, P);
fprintf('reference sm SD (/10) = %.3f\n', sd_ref);

% ── loop over drive conditions: build drive -> generate -> analyse ──────────
nC = numel(conditions);
R = struct('name',{},'mc',{},'res',{},'info',{},'acaus_ratio',{});
fprintf('\nComparing drive constructions ...\n');
for c = 1:nC
    name = conditions{c};
    rng(P.seed);                                   % same synthetic-noise seed across runs
    [mc, sm_in] = build_drive(name, motion_cache, bl_t, P, sd_ref, tau_ac, fps);

    Pm = P;
    if calibrate, Pm.theta = calibrate_theta(bl_t, sm_in, Pm, real_med); end

    [bl, info] = generate_accumulator_dataset(bl_t, sm_in, Pm);

    opts.fig_title = name;
    res = run_hazard_kernel(bl, mc, opts);

    cm = res.caus_mask;
    cmag = mean(abs(res.kernel(cm)));              % mean |beta| over causal lags
    amag = mean(abs(res.kernel(~cm)));             % ... over acausal (future) lags
    ratio = amag / max(cmag, eps);

    R(c).name = name; R(c).mc = []; R(c).res = res; R(c).info = info;
    R(c).acaus_ratio = ratio;

    fprintf(['%-11s | theta=%5.2f | med=%4.0f cens=%4.1f%% | events=%d | ' ...
             '|beta| causal=%.4f acausal=%.4f  (acausal/causal=%.2f)\n'], ...
        name, Pm.theta, info.median_fr, 100*info.cens_frac, res.n_events, cmag, amag, ratio);
end

% ── Figure 1: recovered kernels (RAW, native units), acausal region shaded ──
% Plotting raw beta (NOT shape-normalised): a per-curve max-norm makes a small,
% flat kernel look as tall as a large causal one. Native units show the true
% amplitude difference (level-only kernels are genuinely small).
cmap = lines(nC);
kmax = max(arrayfun(@(c) max(abs(R(c).res.kernel)), 1:nC));   % data-driven y-limit
yl   = 1.15 * kmax * [-1 1];
figure('Color','w','Position',[80 80 860 480]); hold on
patch([0 opts.kernel_future_s opts.kernel_future_s 0], [yl(1) yl(1) yl(2) yl(2)], ...
      [0.92 0.92 0.92], 'EdgeColor','none');          % acausal control region (tau>0)
text(opts.kernel_future_s/2, 0.95*yl(2), 'acausal (future)', 'HorizontalAlignment','center', 'FontSize',9)
hh = gobjects(1,nC);
for c = 1:nC
    hh(c) = plot(R(c).res.t_rel, R(c).res.kernel, 'Color',cmap(c,:), 'LineWidth',2);
end
yline(0,'k:'); xline(0,'--k','escape'); ylim(yl)
xlabel('Time relative to un-freeze (s)   —   past < 0 < future')
ylabel('Recovered kernel  \beta(\tau)')
title('Kernel by drive construction (raw \beta, native units; flat & symmetric = level-only)')
legend(hh, conditions, 'Location','northwest', 'box','off')

% ── Figure 2: acausal/causal weight ratio (the discriminator) ───────────────
figure('Color','w','Position',[80 80 480 420]); hold on
bar(1:nC, [R.acaus_ratio], 0.6, 'FaceColor',[0.4 0.6 0.85], 'EdgeColor','none')
yline(1, 'r--', 'level-only (\approx symmetric)')
set(gca,'XTick',1:nC,'XTickLabel',conditions)
ylabel('acausal / causal  mean |\beta(\tau)|')
title('Temporal-structure test: low = causal/time-resolved, \approx1 = level-only')

% ── Figure 3 (optional): raw kernel vs event-triggered sm excess (eta_exc) ──
% Tells whether kernel features track REAL event-triggered sm structure (e.g. a
% post-escape rise in the acausal/future region) rather than fit/edge artifacts:
% if the kernel and eta_exc rise together at the same lags, the feature is real.
if do_eta_overlay
    figure('Color','w','Position',[80 80 460*nC 420]);
    for c = 1:nC
        t = R(c).res.t_rel; k = R(c).res.kernel; ee = R(c).res.eta_exc;
        subplot(1,nC,c); hold on
        kl = 1.15*max(abs(k))*[-1 1];
        patch([0 opts.kernel_future_s opts.kernel_future_s 0], [kl(1) kl(1) kl(2) kl(2)], ...
              [0.92 0.92 0.92], 'EdgeColor','none');
        yyaxis left
        plot(t, k, 'LineWidth',2); ylabel('kernel  \beta(\tau)'); ylim(kl)
        yyaxis right
        plot(t, ee, '--', 'LineWidth',1.4); ylabel('event-triggered sm excess  \eta_{exc}(\tau)')
        xline(0,'--k','escape')
        xlabel('Time rel. to un-freeze (s)')
        title(R(c).name)
    end
    sgtitle('Raw kernel (left axis) vs event-triggered sm excess (right axis, dashed)')
end

% ── save ─────────────────────────────────────────────────────────────────────
if save_out
    stamp  = datestr(now, 'yyyymmdd_HHMMSS'); %#ok<TNOW1,DATST>
    outdir = fullfile(here, 'results');
    if ~exist(outdir, 'dir'), mkdir(outdir); end
    save(fullfile(outdir, ['kernel_temporal_structure_' stamp '.mat']), ...
        'R', 'conditions', 'P', 'tau_gen', 'tau_ac', 'opts', 'real_med', 'sd_ref', '-v7.3');
    fprintf('\nSaved temporal-structure comparison to %s (stamp %s)\n', outdir, stamp);
end

% ════════════════════════════════════════════════════════════════════════
% Local functions
% ════════════════════════════════════════════════════════════════════════
function [mc, sm_in] = build_drive(name, motion_cache, bl_t, P, sd_ref, tau_ac, fps)
% Return the motion_cache the ANALYSIS reads (mc) and the per-bout generative
% drive (sm_in), consistently, for one drive construction.
    nB = height(bl_t);
    switch name
        case {'timeseries','mean'}
            mc    = motion_cache;                     % analysis reads the REAL sm
            scale = 1/10;                             % real-sm convention
        case {'white','ar1'}
            mc = containers.Map('KeyType','double','ValueType','any');
            flies = unique(bl_t.fly);
            for i = 1:numel(flies)
                L = numel(motion_cache(flies(i)));
                switch name
                    case 'white'
                        x = sd_ref * randn(L,1);
                    case 'ar1'
                        phi = exp(-1/(max(tau_ac,eps)*fps));
                        x   = filter(sd_ref*sqrt(1-phi^2), [1 -phi], randn(L,1));  % stationary SD ~ sd_ref
                end
                mc(double(flies(i))) = x;
            end
            scale = 1;                                % synthetic already at sd_ref
        otherwise
            error('unknown drive condition: %s', name);
    end

    sm_in = cell(nB,1);
    for b = 1:nB
        v = scale * local_win(mc(bl_t.fly(b)), bl_t.onsets(b), P.maxT_fr);
        v(isnan(v)) = 0;
        sm_in{b} = v;
    end

    % 'mean': collapse each bout's drive to its scalar mean (constant drift) —
    % the analysis still sees the time-resolved mc above.
    if strcmp(name, 'mean')
        sm_in = cellfun(@(v) mean(v) * ones(size(v)), sm_in, 'UniformOutput', false);
    end
end

function w = local_win(sm, on, maxT_fr)
    sm  = sm(:);
    cap = min(maxT_fr, numel(sm) - on + 1);
    w   = sm(on : on + cap - 1);
end

function sd = ref_sm_sd(motion_cache, bl_t, P)
    nB = height(bl_t); v = cell(nB,1);
    for b = 1:nB
        v{b} = local_win(motion_cache(bl_t.fly(b))/10, bl_t.onsets(b), P.maxT_fr);
    end
    sd = std(cell2mat(cellfun(@(x) x(:), v, 'UniformOutput', false)), 'omitnan');
end

function theta = calibrate_theta(bl_t, sm_in, P, target_med_fr)
% Bisection on theta so the median synthetic (uncensored, truncated) duration
% ~ the real median (higher theta -> longer durations; monotone).
    lo = 0.02; hi = 200; theta = NaN;
    Pc = P; Pc.n_aug = 1;
    for it = 1:40 %#ok<NASGU>
        theta    = 0.5*(lo + hi);
        Pc.theta = theta;
        [~, info] = generate_accumulator_dataset(bl_t, sm_in, Pc);
        if numel(info.sim_durs) < 50, m = NaN; else, m = info.median_fr; end
        if isnan(m),                hi = theta;
        elseif m < target_med_fr,   lo = theta;
        else,                       hi = theta;
        end
        if (hi - lo) < 1e-3, break; end
    end
end
