clearvars; close all; clc;
% ════════════════════════════════════════════════════════════════════════
% KERNEL RECOVERY SWEEP
%
% Generate synthetic freeze datasets from a leaky accumulator that integrates
% the REAL social-motion each fly experienced, AUGMENTING every real bout into
% P.n_aug synthetic trajectories. Run the hazard-kernel analysis
% (run_hazard_kernel — the function-form core of hazard_kernel_sm.m, IDENTICAL
% methodology) on each synthetic dataset and ask whether it recovers what we
% put in. Two sweep modes (set sweep_var):
%
%   'beta' — social GAIN, leak fixed (specificity / amplitude recovery):
%            beta = 0  => sm does not drive the accumulator => flat/null kernel.
%            Recovered kernel AMPLITUDE should grow with beta and vanish at 0.
%   'tau'  — leak TIME-CONSTANT, gain fixed (timescale recovery):
%            tau = Inf -> flat kernel; tau finite -> exponential kernel ~ tau;
%            tau -> dt -> delta kernel. Recovered tau_com/tau_exp should track.
%
% Because the drive is the SAME real sm the analysis later reads back, the
% regressor autocorrelation matches the real analysis — an honest parameter-
% recovery test, not a toy.
% ════════════════════════════════════════════════════════════════════════

% ── paths ──────────────────────────────────────────────────────────────────
here     = fileparts(mfilename('fullpath'));
code_dir = fileparts(fileparts(here));                       % .../code
addpath(fullfile(code_dir, 'sims', 'simulators'));            % sim_leaky_accumulator
addpath(fullfile(code_dir, 'social_motion', 'proximity'));    % run_hazard_kernel

% ── generative regime (drive = mu0 + beta*sm) ──────────────────────────────
fps           = 60;
P.dt          = 1/fps;
P.mu0         = 3;      % baseline drive  -> sets the ~real-median duration
P.beta        = 1;      % social gain (used as the FIXED gain when sweeping tau)
P.sigma       = 1.0;    % diffusion SD
P.seed        = 1000;
P.maxT_fr     = 630;    % 10.5 s integration cap (no-crossing => censored here)
P.trunc_point = 15;     % min duration (mirror real left-truncation)
P.n_aug       = 2;      % synthetic trajectories per real bout (augmentation)
P.theta       = 10;      % bound (fixed unless calibrate = true)

% ── what to sweep ───────────────────────────────────────────────────────────
sweep_var  = 'beta';            % 'beta' (social gain) or 'tau' (leak time-const)
sweep_list = [5 3 1 0];         % swept values: beta (gain) or tau (seconds)
tau_fixed  = Inf;               % leak used when sweeping beta (Inf = perfect accumulator)

calibrate  = false;   % per-value bisection of theta -> real median duration
save_out   = true;

% ── hazard-analysis options (match hazard_kernel_sm.m methodology) ──────────
% n_shuffle>0 and do_stage2=true add the temporal-order control and the
% integration-vs-extrema model comparison, but are slow — off by default.
opts = struct('fps',fps, 'kernel_past_s', 3, 'kernel_future_s', 2, 'dt_frames', 6, ...
    'nb_kernel', 9, 'nb_acausal', 6, 'bump_width', 2, 'nb_baseline',6, ...
    'cv_folds', 5, 'n_shuffle', 0, 'shuffle_mode', 'full', 'trunc_point', P.trunc_point, ...
    'do_plot',true, 'do_stage2', false, 'use_penalty', true, 'anchor_zero',true);

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

% real median duration (theta calibration target + KS reference)
if ismember('durations', bl_t.Properties.VariableNames)
    real_dur = double(bl_t.durations);
else
    real_dur = double(bl_t.ends - bl_t.onsets);
end
real_med  = median(real_dur);
real_durs = real_dur(real_dur > P.trunc_point & real_dur < P.maxT_fr);
fprintf('Template: %d real freezes, median %.0f fr (%.2f s)\n', height(bl_t), real_med, real_med/fps);

% ── per-bout sm drive: onset forward, capped, kept >= 0 (NOT normalized) ────
% Read motion_cache DIRECTLY (the same series run_hazard_kernel reads back), so
% the generative drive and the analysed regressor are the same signal. The
% accumulator integrates social motion in its RAW units (here only your /10);
% P.beta sets the gain directly. We do NOT z-score / SD-scale the drive: that
% would merely rescale beta, and run_hazard_kernel z-scores sm internally, so it
% cannot change the recovered kernel SHAPE anyway.
nB = height(bl_t); sm_in = cell(nB,1);
for b = 1:nB
    sm  = motion_cache(bl_t.fly(b)); sm = sm(:) ./ 10;
    on  = bl_t.onsets(b);
    cap = min(P.maxT_fr, numel(sm) - on + 1);
    sm_in{b} = sm(on : on + cap - 1);
    sm_in{b}(isnan(sm_in{b})) = 0;        % missing sm -> baseline drive only
end
sm_sd = std(cell2mat(cellfun(@(v) v(:), sm_in, 'UniformOutput', false)), 'omitnan');
fprintf('sm drive: SD = %.3f (raw /10 units, NOT normalized)\n', sm_sd);

% ── sweep: (calibrate theta) -> generate -> analyse ─────────────────────────
R = struct('val',{},'beta',{},'tau',{},'lambda',{},'theta',{},'info',{},'res',{},'ks',{},'bl',{});
fprintf('\nSweeping %s; generating + analysing each value ...\n', sweep_var);
for m = 1:numel(sweep_list)
    v  = sweep_list(m);
    Pm = P;
    switch sweep_var
        case 'beta', Pm.beta = v;  tau = tau_fixed;
        case 'tau',                tau = v;
        otherwise, error('sweep_var must be ''beta'' or ''tau''');
    end
    Pm.lambda = 0; if ~isinf(tau), Pm.lambda = 1/tau; end

    if calibrate
        Pm.theta = calibrate_theta(bl_t, sm_in, Pm, real_med);
    end

    [bl, info] = generate_accumulator_dataset(bl_t, sm_in, Pm);

    ks = NaN;
    if numel(info.sim_durs) > 50, [~,~,ks] = kstest2(info.sim_durs, real_durs); end

    opts.fig_title = sweep_title(sweep_var, v, tau);
    res = run_hazard_kernel(bl, motion_cache, opts);

    R(m).val=v; R(m).beta=Pm.beta; R(m).tau=tau; R(m).lambda=Pm.lambda; R(m).theta=Pm.theta;
    R(m).info=info; R(m).res=res; R(m).ks=ks; R(m).bl=bl;

    fprintf(['%s=%6.3g | theta=%5.2f | med=%4.0f cens=%4.1f%% filt=%4.1f%% KS=%.3f | ' ...
             'tau_com=%.2f tau_exp=%.2f R2=%.2f | events=%d\n'], ...
        sweep_var, v, Pm.theta, info.median_fr, 100*info.cens_frac, 100*info.filt_frac, ks, ...
        res.tau_com, res.tau_hat, res.R2_exp, res.n_events);
end

xval = [R.val];
cmap = parula(numel(R));
lab  = arrayfun(@(k) sweep_label(sweep_var, R(k).val), 1:numel(R), 'UniformOutput', false);

% per-value kernel amplitude (peak and net causal weight)
kpk  = arrayfun(@(k) max(R(k).res.kernel(R(k).res.caus_mask)), 1:numel(R));
ksum = arrayfun(@(k) sum(R(k).res.kernel(R(k).res.caus_mask)), 1:numel(R));

% ── Figure 1: recovered kernels across the sweep ────────────────────────────
figure('Color','w','Position',[80 80 820 480]); hold on
hh = gobjects(1,numel(R));
for m = 1:numel(R)
    rr = R(m).res; cm = rr.caus_mask;
    hh(m) = plot(rr.t_rel(cm), rr.kernel(cm), 'Color',cmap(m,:), 'LineWidth',1.8);
end
yline(0,'k:'); xline(0,'--k','escape');
xlabel('Time relative to un-freeze (s)   (causal: \tau \leq 0)')
ylabel('Recovered kernel \beta(\tau)  (per SD of sm)')
title(sprintf('Recovered hazard kernels across the %s sweep', sweep_var))
legend(hh, lab, 'Location','northwest', 'box','off')

% ── Figure 2: recovery readout (depends on what we swept) ───────────────────
figure('Color','w','Position',[80 80 560 470]); hold on
switch sweep_var
    case 'beta'
        % amplitude recovery: kernel strength should grow with beta, ~0 at beta=0
        yline(0,'k:');
        plot(xval, kpk,  'o-', 'LineWidth',1.6, 'DisplayName','peak causal \beta(\tau)')
        plot(xval, ksum, 's--','LineWidth',1.2, 'DisplayName','\Sigma causal \beta(\tau)')
        xlabel('Generative social gain \beta'); ylabel('Recovered kernel amplitude')
        title('Gain recovery: amplitude vs \beta   (\beta=0 \Rightarrow flat/null)')
        legend('Location','northwest', 'box','off')
    case 'tau'
        taus = xval; finite = ~isinf(taus);
        tcom = arrayfun(@(k) R(k).res.tau_com, 1:numel(R));
        texp = arrayfun(@(k) R(k).res.tau_hat, 1:numel(R));
        if any(finite)
            lim = [min(taus(finite))/2, max(taus(finite))*2];
            plot(lim, lim, 'k:', 'LineWidth',1)                            % identity
            plot(taus(finite), tcom(finite), 'o-', 'LineWidth',1.6, 'DisplayName','\tau_{com}')
            plot(taus(finite), texp(finite), 's--','LineWidth',1.2, 'DisplayName','\tau_{exp}')
        end
        if any(~finite), yline(tcom(~finite), ':', 'accum \tau_{com}', 'Color',[.4 .4 .4]); end
        set(gca,'XScale','log','YScale','log'); grid on
        xlabel('Generative leak \tau (s)'); ylabel('Recovered \tau (s)')
        title('Time-horizon recovery (identity line = perfect)')
        legend('Location','northwest', 'box','off')
end

% ── Figure 3: duration-distribution match (truncated CDFs vs real) ──────────
figure('Color','w','Position',[80 80 700 470]); hold on
[fr,xr] = ecdf(real_durs); stairs(xr/fps, fr, 'k', 'LineWidth',2.5, 'DisplayName','real');
for m = 1:numel(R)
    if numel(R(m).info.sim_durs) > 50
        [fs,xs] = ecdf(R(m).info.sim_durs);
        stairs(xs/fps, fs, 'Color',cmap(m,:), 'LineWidth',1, 'DisplayName',lab{m});
    end
end
xlabel('Freeze duration (s)'); ylabel('CDF'); xlim([0 P.maxT_fr/fps])
title('Truncated duration distributions: real (black) vs each value')
legend('Location','southeast', 'box','off')

% ── save sweep (synthetic datasets + recovery results + params) ─────────────
if save_out
    stamp  = datestr(now, 'yyyymmdd_HHMMSS'); %#ok<TNOW1,DATST>
    outdir = fullfile(here, 'results');
    if ~exist(outdir, 'dir'), mkdir(outdir); end
    save(fullfile(outdir, ['kernel_recovery_' sweep_var '_' stamp '.mat']), ...
        'R', 'sweep_var', 'sweep_list', 'tau_fixed', 'P', 'opts', ...
        'real_med', 'real_durs', 'sm_sd', '-v7.3');
    fprintf('\nSaved %s-sweep (datasets + recovery) to %s (stamp %s)\n', sweep_var, outdir, stamp);
end

% ════════════════════════════════════════════════════════════════════════
% Local functions
% ════════════════════════════════════════════════════════════════════════
function theta = calibrate_theta(bl_t, sm_in, P, target_med_fr)
% Bisection on theta so the median synthetic (uncensored, truncated) duration
% ~ the real median. Higher theta -> harder to cross -> longer durations
% (monotone), so the bracket logic is well-posed. One augmentation rep is
% enough to estimate the median.
    lo = 0.02; hi = 200; theta = NaN;
    Pc = P; Pc.n_aug = 1;
    for it = 1:40 %#ok<NASGU>
        theta   = 0.5*(lo + hi);
        Pc.theta = theta;
        [~, info] = generate_accumulator_dataset(bl_t, sm_in, Pc);
        if numel(info.sim_durs) < 50, m = NaN; else, m = info.median_fr; end
        if isnan(m)
            hi = theta;                 % too few crossings -> theta too high
        elseif m < target_med_fr
            lo = theta;                 % too fast -> raise theta
        else
            hi = theta;                 % too slow -> lower theta
        end
        if (hi - lo) < 1e-3, break; end
    end
end

function s = sweep_label(sweep_var, v)
    switch sweep_var
        case 'beta', s = sprintf('\\beta=%.2g', v);
        case 'tau'
            if isinf(v), s = '\tau=\infty (accum)'; else, s = sprintf('\\tau=%.2g s', v); end
    end
end

function s = sweep_title(sweep_var, v, tau)
    switch sweep_var
        case 'beta'
            if isinf(tau), ts = '\tau=\infty'; else, ts = sprintf('\\tau=%.2g s', tau); end
            s = sprintf('\\beta=%.2g  (%s)', v, ts);
        case 'tau'
            if isinf(v), s = 'accumulator (\tau=\infty)'; else, s = sprintf('\\tau=%.2g s', v); end
    end
end
