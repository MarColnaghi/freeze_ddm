clearvars; clc; close all;
% ════════════════════════════════════════════════════════════════════════
% MEMORY-TIMESCALE ANALYSIS — Steps 1–5 of the hierarchical GLM pipeline.
%
% Inverts the inference: instead of recovering a temporal kernel and reading a
% leak off its shape (multicollinear, ringy, weakly identified), we ask which
% single causal summary of recent social motion best predicts un-freezing, and
% over what timescale. Each model is ONE z-scored predictor on top of a shared
% baseline hazard + controls, scored by bout-grouped cross-validated logLik.
%
%   Step 1  motion matters         : baseline      vs  baseline + mean motion
%   Step 2  temporal structure     : mean          vs  mean + slope (derivative)
%   Step 3  boxcar timescale       : CV-LL( window T )      → averaging horizon
%   Step 4  exponential timescale  : CV-LL( leak τ )        → leak time-constant
%   Step 5  memory architecture    : instantaneous vs best boxcar vs best exp
%
% Steps 6–9 (kernel confirmation, bootstrap, shuffle nulls, DDM) are separate
% scripts that reuse this engine and the SAME design/CV folds.
% ════════════════════════════════════════════════════════════════════════

addpath('/Users/marcocolnaghi/PhD/freeze_ddm/code/social_motion/proximity');   % path_generator, data parsing
here = fileparts(mfilename('fullpath')); addpath(here);                          % glm engine

% ── Settings ──────────────────────────────────────────────────────────────
dt_frames   = 6;        % hazard sampled every 6 frames (0.1 s)
nb_baseline = 6;
cv_folds    = 5;
trunc_point = 30;       % 0.5 s min duration (left truncation / delayed entry)
contact_threshold = 80; % px — right-censor a freeze at first contact-distance frame
slope_win_s = 1.0;      % window for the Step-2 derivative feature
T_grid   = [0.25 0.5 1 2 4 8];   % boxcar windows (s)
tau_grid = [0.25 0.5 1 2 4 8];   % exponential leak timescales (s)
save_out = true;
rng(7);

dopts = struct('fps',60,'dt_frames',dt_frames,'nb_baseline',nb_baseline, ...
    'cv_folds',cv_folds,'trunc_point',trunc_point,'rng_seed',7);
fps = dopts.fps;

% ── Load + filter (same pipeline as hazard_kernel_sm.m) ───────────────────
threshold_imm = 2; threshold_mob = 2; threshold_pc = 4;
id_code = sprintf('imm%d_mob%d_pc%d', threshold_imm, threshold_mob, threshold_pc);
paths        = path_generator('folder','social_motion/proximity','bouts_id',id_code,'imfirst',false);
bouts        = importdata(fullfile(paths.dataset,'bouts.mat'));
motion_cache = importdata(fullfile(paths.cache_path,'motion_cache.mat'));
thresholds   = define_thresholds('le_window', struct('le_window_sl',[5 55],'le_window_fl',[5 55]));
bouts        = bouts_formatting(bouts, thresholds);

bl = data_parser_new(bouts,'type','immobility','period','loom','window','le', ...
        'nloom',2:20,'min_dur',trunc_point);
bl = impose_contact_threshold(bl, 'threshold', contact_threshold, 'type', 'onlyfreeze');
bl.dur_frames = bl.ending_time;
bl.bout_id    = (1:height(bl))';
bl = bl(bl.dur_frames > trunc_point, :);

% ── Shared design (rows / events / baseline / controls / CV folds) ────────
D = build_hazard_design(bl, motion_cache, dopts);

R = struct();   % collected results

% ════════════════════════════════════════════════════════════════════════
% STEP 1 — does social motion matter at all?
%   baseline  vs  baseline + mean motion (1 s boxcar as the canonical "level")
% ════════════════════════════════════════════════════════════════════════
fprintf('\n══ STEP 1: motion matters ══\n');
m_base = fit_cv(D, [], {});
f_mean = memory_feature(D, motion_cache, struct('type','box','timescale',1.0));
m_mean = fit_cv(D, f_mean, {'mean'});
R.step1.dCVLL = m_mean.cvll - m_base.cvll;
R.step1.dAIC  = m_base.aic  - m_mean.aic;     % positive = motion model better
fprintf('baseline       CV-LL=%.5f  AIC=%.1f\n', m_base.cvll, m_base.aic);
fprintf('+ mean motion  CV-LL=%.5f  AIC=%.1f\n', m_mean.cvll, m_mean.aic);
fprintf('ΔCV-LL = %+.5f   ΔAIC = %+.1f  (>10 ⇒ motion predicts unfreezing)\n', ...
    R.step1.dCVLL, R.step1.dAIC);

% ════════════════════════════════════════════════════════════════════════
% STEP 2 — does temporal STRUCTURE (the derivative) matter beyond level?
%   mean  vs  mean + slope
% ════════════════════════════════════════════════════════════════════════
fprintf('\n══ STEP 2: temporal structure / derivative ══\n');
f_slope = memory_feature(D, motion_cache, struct('type','slope','timescale',slope_win_s));
ok = ~isnan(f_slope);   % slope is NaN only for early-recording edges (rare)
m_ms = fit_cv(D, [f_mean, f_slope], {'mean','slope'});
R.step2.dCVLL = m_ms.cvll - m_mean.cvll;
R.step2.dAIC  = m_mean.aic - m_ms.aic;
fprintf('mean           CV-LL=%.5f  AIC=%.1f\n', m_mean.cvll, m_mean.aic);
fprintf('mean + slope   CV-LL=%.5f  AIC=%.1f\n', m_ms.cvll, m_ms.aic);
fprintf('ΔCV-LL = %+.5f   ΔAIC = %+.1f  (slope coeff p=%.3g; %d/%d rows have slope)\n', ...
    R.step2.dCVLL, R.step2.dAIC, m_ms.mdl.Coefficients.pValue(strcmp(m_ms.mdl.CoefficientNames,'slope')), ...
    sum(ok), numel(ok));

% ════════════════════════════════════════════════════════════════════════
% STEP 3 — model-free averaging timescale: CV-LL vs boxcar window T
% ════════════════════════════════════════════════════════════════════════
fprintf('\n══ STEP 3: boxcar memory timescale ══\n');
Pbox = profile_timescale(D, motion_cache, 'box', T_grid);
R.box = Pbox;
for i = 1:numel(T_grid)
    fprintf('  T=%5.2f s : CV-LL=%.5f ± %.5f\n', T_grid(i), Pbox.cvll(i), Pbox.cvse(i));
end
fprintf('best boxcar window ≈ %.2f s\n', Pbox.best_ts);

% ════════════════════════════════════════════════════════════════════════
% STEP 4 — leak timescale: CV-LL vs exponential (leaky-integrator) τ
% ════════════════════════════════════════════════════════════════════════
fprintf('\n══ STEP 4: exponential (leak) timescale ══\n');
Pexp = profile_timescale(D, motion_cache, 'exp', tau_grid);
R.exp = Pexp;
for i = 1:numel(tau_grid)
    fprintf('  τ=%5.2f s : CV-LL=%.5f ± %.5f\n', tau_grid(i), Pexp.cvll(i), Pexp.cvse(i));
end
fprintf('best leak τ ≈ %.2f s\n', Pexp.best_ts);

% ════════════════════════════════════════════════════════════════════════
% STEP 5 — memory architecture: instantaneous vs best boxcar vs best exp
%   (matched complexity: one predictor each; CV-LL is scale-free)
% ════════════════════════════════════════════════════════════════════════
fprintf('\n══ STEP 5: memory architecture ══\n');
f_inst = memory_feature(D, motion_cache, struct('type','inst','timescale',0));
m_inst = fit_cv(D, f_inst, {'mem'});
f_box  = memory_feature(D, motion_cache, struct('type','box','timescale',Pbox.best_ts));
m_box  = fit_cv(D, f_box, {'mem'});
f_exp  = memory_feature(D, motion_cache, struct('type','exp','timescale',Pexp.best_ts));
m_exp  = fit_cv(D, f_exp, {'mem'});
arch = struct('name', {'instantaneous', sprintf('boxcar %.2gs',Pbox.best_ts), sprintf('exp τ=%.2gs',Pexp.best_ts)}, ...
              'cvll', {m_inst.cvll, m_box.cvll, m_exp.cvll}, ...
              'aic',  {m_inst.aic,  m_box.aic,  m_exp.aic});
R.arch = arch;
fprintf('%-18s %12s %10s\n','architecture','CV-LL','AIC');
for i = 1:numel(arch)
    fprintf('%-18s %12.5f %10.1f\n', arch(i).name, arch(i).cvll, arch(i).aic);
end
[~,bi] = max([arch.cvll]);
fprintf('best architecture by CV-LL: %s\n', arch(bi).name);
fprintf('exp − box ΔCV-LL = %+.5f  (>0 ⇒ leak beats finite window)\n', m_exp.cvll - m_box.cvll);

% ════════════════════════════════════════════════════════════════════════
% Figure
% ════════════════════════════════════════════════════════════════════════
fh = figure('Color','w','Position',[60 60 1180 760]);
tiledlayout(2,2,'TileSpacing','compact','Padding','compact')

nexttile;   % Steps 1–2 bar of ΔCV-LL
bar([R.step1.dCVLL, R.step2.dCVLL]); grid on
set(gca,'XTickLabel',{'+mean vs base','+slope vs mean'})
ylabel('\Delta CV-LL / obs'); title('Steps 1–2: motion & its derivative add predictive power')

nexttile;   % Step 3 boxcar profile
errorbar(T_grid, Pbox.cvll, Pbox.cvse, 'o-','LineWidth',1.5); set(gca,'XScale','log'); grid on; hold on
plot(Pbox.best_ts, Pbox.best_cvll, 'rp','MarkerFaceColor','r','MarkerSize',12)
xlabel('boxcar window T (s)'); ylabel('CV-LL / obs'); title('Step 3: averaging timescale (boxcar)')

nexttile;   % Step 4 exponential profile
errorbar(tau_grid, Pexp.cvll, Pexp.cvse, 'o-','LineWidth',1.5,'Color',[0.85 0.4 0.2]); set(gca,'XScale','log'); grid on; hold on
plot(Pexp.best_ts, Pexp.best_cvll, 'rp','MarkerFaceColor','r','MarkerSize',12)
xlabel('leak \tau (s)'); ylabel('CV-LL / obs'); title('Step 4: leak timescale (exponential)')

nexttile;   % Step 5 architecture bars (relative to baseline)
vals = [arch.cvll] - m_base.cvll;
bar(vals); grid on; set(gca,'XTickLabel',{arch.name})
ylabel('CV-LL − baseline'); title('Step 5: memory architecture')

% ════════════════════════════════════════════════════════════════════════
% Save
% ════════════════════════════════════════════════════════════════════════
if save_out
    stamp = datestr(now,'yyyymmdd_HHMMSS'); %#ok<TNOW1,DATST>
    params = struct('dt_frames',dt_frames,'nb_baseline',nb_baseline,'cv_folds',cv_folds, ...
        'trunc_point',trunc_point,'contact_threshold',contact_threshold, ...
        'slope_win_s',slope_win_s,'T_grid',T_grid,'tau_grid',tau_grid);
    try
        save(fullfile(paths.dataset, ['glm_memory_timescale_' stamp '.mat']), ...
            'R','params','T_grid','tau_grid','-v7.3');
        exportgraphics(fh, fullfile(paths.dataset, ['glm_memory_timescale_' stamp '.png']), 'Resolution',200);
        fprintf('\nSaved Steps 1–5 (stamp %s)\n', stamp);
    catch ME
        warning('save failed: %s', ME.message);
    end
end
