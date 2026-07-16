
%% ════════════════════════════════════════════════════════════════════════
%  Simulator consistency checks
%  --------------------------------------------------------------------------
%  (1) leaky(lambda=0)     ==  drift_diff_new          trial-by-trial (same seed
%      -> same randn draws -> identical first-passage frame / RT).
%  (2) leaky(lambda)       ~=  sim_ddm_seeded          distributionally (different
%      RNGs -> NOT trial-by-trial identical, but the SAME FPT distribution).
%  (3) trajectory overlap  (same seed) at both ends of the leak continuum.
%  (4) leaky(lambda=1/dt)  ==  extrema_detection_clean trial-by-trial.
%  (5) leak sweep: PERFECT integrator (lambda=0, flat memory) -> EXTREMA
%      detector (lambda=1/dt, delta memory), parametrically.
%  Shared dynamics: x(k+1)=x(k)(1-lambda*dt)+drift(k)*dt+sigma*sqrt(dt)*randn,
%  start x=0, cross at x>=theta, drift = per-frame9 RATE, frame k <-> time k*dt.
% ════════════════════════════════════════════════════════════════════════
this_dir = fileparts(mfilename('fullpath'));
addpath(fullfile(this_dir, 'simulators'));

dt    = 1/60;                          % s / frame (native social-motion rate)
sigma = 1.0;                           % diffusion SD
theta = 0.35;                           % bound -- keep SMALL: at the extrema end
                                       % (checks 4/5) a crossing must happen
                                       % within ONE frame, where the noise SD is
                                       % only sigma*sqrt(dt) ~ 0.13, so a large
                                       % theta would make the extrema detector
                                       % ~never fire (endpoint -> 100% censored).
beta  = 3.12;                           % drift gain (rate per unit social motion)

% ── Drive: a REAL per-frame social-motion signal (one bout), loaded with the
%    same pipeline as the validation section at the bottom of this file.
%    Set use_real=false to fall back to a synthetic time-varying drift.
use_real = true;
if use_real
    paths      = path_generator('folder', 'sims/sanity_checks');
    idx_trial  = 950; 2021; 5222; 11;3394; 2020;1918;% 1001;% 100 is great; 4548; 878; 917; 6666;3945;412;%7327;%3331;%2305;%2132; %3331 is fantastic;%3945;
    bouts      = importdata(fullfile(paths.dataset, 'bouts.mat'));
    bouts_proc = data_parser_new(bouts, 'type','immobility', 'period','loom', ...
                                 'window','le', 'nloom', 2:20);
    bouts_proc.smp       = ones(height(bouts_proc), 1);
    bouts_proc.intercept = ones(height(bouts_proc), 1);
    sm_signal  = extract_sm_from_bouts(bouts_proc, 'type','onsets', ...
                                 'output_type','mat', 'window', [0 630], 'norm_factor', 10);
    drift = sm_signal(idx_trial, :).' * beta;   % per-frame drift RATE (column)
    drift = fillmissing(drift, 'previous');      % guard against motion gaps
else
    tt    = (1:600)' * dt;
    drift = 0.5 + 0.4*sin(2*pi*0.2*tt);          % synthetic time-varying rate
end
n     = numel(drift);                  % frames in the drive
T_dd  = n * dt;                        % so drift_diff_new's 0:dt:T grid has n steps

%% ── (1) IDENTITY: leaky(lambda=0) vs drift_diff_new (same seed) ──────────
% NB: exact identity hinges on drift_diff_new drawing its noise as ONE
%     randn(n,1) batch -- the same sequence the leaky loop draws frame by
%     frame. If that draw ever changes, this quietly degrades to a
%     distributional (not trial-by-trial) match.
n_id     = 2000;
frame_la = nan(n_id,1);  rt_la = nan(n_id,1);
frame_dd = nan(n_id,1);  rt_dd = nan(n_id,1);
for s = 1:n_id

    k = sim_leaky_accumulator(drift, theta, sigma, 0, dt, s);
    frame_la(s) = k;  rt_la(s) = k * dt;

    r = drift_diff_new('mu_t', drift, 'theta', theta, 'sigma', sigma, ...
                       'dt', dt, 'z', 0, 'T', T_dd, 'ndt', 0, 'seed', s);
    rt_dd(s) = r;
    if ~isnan(r), frame_dd(s) = round(r / dt); end   % grid index-1 -> frame count
end

both_cross     = ~isnan(frame_la) & ~isnan(frame_dd);
frame_mismatch = sum(both_cross & (frame_la ~= frame_dd));
nan_mismatch   = sum(xor(isnan(frame_la), isnan(frame_dd)));
max_rt_diff    = max(abs(rt_la(both_cross) - rt_dd(both_cross)));
if isempty(max_rt_diff), max_rt_diff = 0; end

fprintf('\n[1] leaky(lambda=0) vs drift_diff_new  (%d seeds)\n', n_id);
fprintf('    crossings (both)          : %d\n', sum(both_cross));
fprintf('    crossing-frame mismatches : %d\n', frame_mismatch);
fprintf('    NaN/cross disagreements   : %d\n', nan_mismatch);
fprintf('    max |RT diff| (s)         : %.3e\n', max_rt_diff);
if frame_mismatch==0 && nan_mismatch==0 && max_rt_diff < 1e-9
    fprintf('    => PASS: trial-by-trial identical.\n');
else
    fprintf('    => FAIL: not identical (see above).\n');
end

%% ── (2) DISTRIBUTION: leaky(lambda=0) vs sim_ddm_seeded ──────────────────
N = 500000;

lambda = 5;

% leaky: N independent trials from one global stream (seed=[] -> no reseed)
rng(1);
rt_leaky = nan(N,1);
tic 
for i = 1:N
    k = sim_leaky_accumulator(drift, theta, sigma, lambda, dt, []);
    if ~isnan(k), rt_leaky(i) = k*dt; end
end
toc

% C++ mex: params = [dt, sigma^2, x0, theta, lambda]; one call -> N trials
params_mex = [dt, sigma^2, 0, theta, lambda];
tic
rt_mex = sim_ddm_seeded(drift, params_mex, N, 1);   % seconds, NaN if no cross
toc

summ  = @(x) [median(x(~isnan(x))), mean(x(~isnan(x))), std(x(~isnan(x))), mean(isnan(x))];
S     = [summ(rt_leaky); summ(rt_mex)];
names = {sprintf('leaky (lam=%.3g)', lambda), 'sim_ddm_seeded'};
fprintf('\n[2] First-passage distribution  (N=%d, lambda=%.3g)\n', N, lambda);
fprintf('    %-18s %8s %8s %8s %8s\n','model','median','mean','std','cens');
for r = 1:2
    fprintf('    %-18s %8.4f %8.4f %8.4f %8.3f\n', names{r}, S(r,1),S(r,2),S(r,3),S(r,4));
end

a = rt_leaky(~isnan(rt_leaky));  b = rt_mex(~isnan(rt_mex));
try
    [h,p,ks] = kstest2(a,b);
    fprintf('    KS leaky vs mex: D=%.4f, p=%.3g  (h=%d)\n', ks, p, h);
catch
    fprintf('    (kstest2 unavailable; skipping KS test)\n');
end

edges = 0-dt/2:dt:T_dd;
figure('Color','w','Position',[10 10 700 300]); hold on
histogram(a, edges, 'Normalization','pdf','EdgeColor','none','FaceAlpha',0.35, ...
    'DisplayName', sprintf('leaky (\\lambda=%.3g)', lambda));
histogram(b, edges, 'Normalization','pdf','EdgeColor','none','FaceAlpha',0.35, ...
    'DisplayName','sim\_ddm\_seeded');
xlabel('first-passage time (s)'); ylabel('pdf'); legend('box','off');
title(sprintf('FPT distributions (\\lambda=%.3g): median leaky=%.3f s, mex=%.3f s', lambda, S(1,1), S(2,1)));

%% ── (3) Trajectory inspection: same seed -> overlapping paths ────────────
% Needs the setup cell at the top (dt, theta, sigma, drift, T_dd).
s = 7;
figure('Color','w','Position',[10 10 950 320]);

% (a) integrator end: leaky(lambda=0) vs drift_diff_new
[k_la, traj_la, t_la] = sim_leaky_accumulator(drift, theta, sigma, 0, dt, s);
[r_dd, traj_dd, t_dd] = drift_diff_new('mu_t', drift, 'theta', theta, 'sigma', sigma, ...
                                       'dt', dt, 'z', 0, 'T', T_dd, 'ndt', 0, 'seed', s);
subplot(1,2,1); hold on
plot(t_dd, traj_dd, 'LineWidth', 1.6, 'DisplayName','drift\_diff\_new');
plot(t_la, traj_la, '--', 'LineWidth', 1.2, 'DisplayName','leaky (\lambda=0)');
yline(theta,'k:','\theta'); xlabel('time (s)'); ylabel('accumulator x');
legend('box','off','Location','northwest');
title(sprintf('Integrator end: cross leaky=%.3f s, dd=%.3f s', k_la*dt, r_dd));

% (b) extrema end: leaky(lambda=1/dt) vs extrema_detection_clean(drift*dt)
[k_le, traj_le, t_le] = sim_leaky_accumulator(drift, theta, sigma, 1/dt, dt, s);
[k_ed, traj_ed, t_ed] = extrema_detection_clean(drift*dt, theta, sigma, dt, s);
subplot(1,2,2); hold on
plot(t_ed, traj_ed, 'LineWidth', 1.6, 'DisplayName','extrema\_detection\_clean');
plot(t_le, traj_le, '--', 'LineWidth', 1.2, 'DisplayName','leaky (\lambda=1/dt)');
yline(theta,'k:','\theta'); xlabel('time (s)'); ylabel('momentary value');
legend('box','off','Location','northwest');
title(sprintf('Extrema end: cross leaky=%g, ed=%g frames', k_le, k_ed));

%% ── (4) IDENTITY at the extrema end: leaky(lambda=1/dt) vs extrema detector
% decay = 1 - lambda*dt -> 0, so the state forgets the past entirely and
% depends only on the current frame: x(k) = drift(k)*dt + sigma*sqrt(dt)*randn.
% extrema_detection_clean takes a LEVEL signal, so feed signal = drift*dt
% (rate -> level) and the two coincide draw-for-draw under a shared seed.

lambda_ext = 1/dt;                 % decay = 1 - (1/dt)*dt = 0
signal_ext = drift * dt;           % per-frame rate -> per-frame level
n_id2 = 2000;
fr_le = nan(n_id2,1);  fr_ed = nan(n_id2,1);
for s = 1:n_id2
    fr_le(s) = sim_leaky_accumulator(drift,       theta, sigma, lambda_ext, dt, s);
    fr_ed(s) = extrema_detection_clean(signal_ext, theta, sigma, dt,         s);
end
bc = ~isnan(fr_le) & ~isnan(fr_ed);
mm = sum(bc & (fr_le ~= fr_ed));
nm = sum(xor(isnan(fr_le), isnan(fr_ed)));
fprintf('\n[4] leaky(lambda=1/dt) vs extrema_detection_clean  (%d seeds)\n', n_id2);
fprintf('    crossings (both)          : %d\n', sum(bc));
fprintf('    crossing-frame mismatches : %d\n', mm);
fprintf('    NaN/cross disagreements   : %d\n', nm);
if mm==0 && nm==0
    fprintf('    => PASS: leak=1/dt reduces the accumulator to the extrema detector.\n');
else
    fprintf('    => FAIL (see above).\n');
end

%% ── (5) CONTINUUM: leak sweeps perfect integrator -> extrema detector ────
% Effective weighting of past drift is geometric: w(lag) = (1 - lambda*dt).^lag.
%   lambda = 0      -> decay=1 -> FLAT  kernel (perfect integrator == drift_diff_new)
%   lambda = 1/dt   -> decay=0 -> DELTA kernel (extrema detector  == extrema_detection_clean)
lam_list = [0 0.5 1 2 5 10 20 40 1/dt];   % 1/dt = 60 here
Nsw      = 200000;

% Shared FPT histogram grid: bin centres at 0, dt, ..., T_dd (one pdf per leak)
edges_sw = -dt/2 : dt : T_dd + dt/2;
ctrs_sw  = 0 : dt : T_dd;
pdf_mat  = nan(numel(ctrs_sw), numel(lam_list));   % time x lambda
med  = nan(numel(lam_list),1);
cens = nan(numel(lam_list),1);
for li = 1:numel(lam_list)
    % distributional only -> use the fast mex (identical dynamics to the .m loop)
    rr  = sim_ddm_seeded(drift, [dt, sigma^2, 0, theta, lam_list(li)], Nsw, 100);
    fin = ~isnan(rr);
    med(li)       = median(rr(fin));
    cens(li)      = mean(~fin);
    pdf_mat(:,li) = histcounts(rr(fin), edges_sw, 'Normalization','pdf')';
end

% Endpoint references from the DEDICATED simulators (smaller N is plenty here):
N_ref = 50000;
rng(100); re = nan(N_ref,1);
for i = 1:N_ref                                                  % extrema reference
    k = extrema_detection_clean(drift*dt, theta, sigma, dt, []);
    if ~isnan(k), re(i) = k*dt; end
end
rd     = sim_ddm_seeded(drift, [dt, sigma^2, 0, theta, 0], N_ref, 7);  % integrator ref
med_dd = median(rd(~isnan(rd)));   % perfect-integrator reference (== drift_diff_new, check 1)
med_ed = median(re(~isnan(re)));   % extrema reference

figure('Color','w','Position',[10 10 900 750]);
tiledlayout(2,2,"TileSpacing", "compact", "Padding", "compact")
color_lambda = cbrewer2('RdYlBu',length(lam_list));

nexttile
hold on            % (a) effective memory kernel vs leak
lags = 0:round(2/dt);              % up to 2 s of lag
for li = 1:numel(lam_list)
    plot(lags*dt, (1 - lam_list(li)*dt).^lags, 'DisplayName', sprintf('\\lambda=%.3g', lam_list(li)), 'LineWidth', 2, 'Color', color_lambda(li,:));
end
legend('Location', 'northoutside', 'FontSize', 12, 'Box', 'off', 'NumColumns', 3)
xlabel('lag (s)'); ylabel('weight on past drift'); ylim([-0.05 1.05]);
apply_generic(gca, 'xlim', [-0.1 2.1], 'font_size', 22)

nexttile
hold on
yyaxis left
plot(lam_list, med, '-','LineWidth', 1.8); ylabel('median RT (s)');
scatter(lam_list, med, 90, color_lambda, 'filled')
ylim([0 5])
yline(med_dd, 'Label','Perfect (\lambda=0)', 'LineStyle',' --', 'FontSize', 12, 'Color',  color_lambda(1,:)); yline(med_ed, 'Label', 'Extrema Detection', 'LineStyle',' --', 'FontSize', 12, 'Color',  color_lambda(end,:));
yyaxis right
ylim([-5 105])
plot(lam_list, 100*cens, '-','LineWidth', 2); ylabel('Censored (%)');
scatter(lam_list, 100*cens, 90, color_lambda, 'filled')

xlabel('\lambda');
apply_generic(gca, 'xlim', [-1.5 1/dt+1.5], 'font_size', 22)


nexttile
imagesc(ctrs_sw, 1:numel(lam_list), pdf_mat');
clim([0 2])
set(gca,'YTick',1:numel(lam_list),'YTickLabel',compose('%.3g',lam_list(:)));
xlim([0 min(T_dd)]);
colormap(cbrewer2('Reds', []))
xlabel('first-passage time (s)'); ylabel('\lambda (1/s)');
apply_generic(gca, 'xlim', [-0.5 10.5], 'font_size', 22)

nexttile
hold on
for li = 1:numel(lam_list)
    plot(ctrs_sw, pdf_mat(:,li), 'Color', [color_lambda(li,:) 0.7], 'LineWidth', 1.2, ...
        'DisplayName', sprintf('\\lambda=%.3g', lam_list(li)));
end
xlim([0 min(T_dd)]); xlabel('first-passage time (s)'); ylabel('pdf');
apply_generic(gca, 'xlim', [-0.5 10.5],'ylim', [0 3], 'font_size', 22)


fprintf('\n[5] Continuum endpoints (median FPT):\n');
fprintf('    leaky(lambda=0)    = %.4f   vs  drift_diff_new          = %.4f\n', med(1),   med_dd);
fprintf('    leaky(lambda=1/dt) = %.4f   vs  extrema_detection_clean = %.4f\n', med(end), med_ed);

%% ─── Original validation: drift_diff_new vs C++ mex on real bout data ────
paths = path_generator('folder', 'sims/sanity_checks');


% Select the parameters of the model
true_bound = 2.9;
beta_drift = 6.2;
true_ndt = 0.5;

% Select the dt for signal upsampling
dt = 1/300;
fps = 1/60;
inizio = 0;
fine = 630;
% Run this now for the timevarying model
points.truncation = 0.5;
points.censoring = 10.5;
n_T = points.censoring / dt;
time_vector = 0:fps:points.censoring;
time_vector_us = 0:dt:points.censoring;

bouts = importdata(fullfile(paths.dataset, 'bouts.mat'));
bouts_proc = data_parser_new(bouts, 'type', 'immobility', 'period', 'loom', 'window', 'le', 'nloom', 2:20);
bouts_proc.smp = ones(height(bouts_proc), 1);
bouts_proc.intercept = ones(height(bouts_proc), 1);

sm = importdata(fullfile(paths.cache_path, "motion_cache.mat"));
sm_signal = extract_sm_from_bouts(bouts_proc, 'type', 'onsets', 'output_type', 'mat', 'window', [inizio fine], 'norm_factor', 10);

sm_signal_us = interp1(time_vector(:), sm_signal', time_vector_us(:), 'nearest');

n_sims = 50000;
mu_tv = sm_signal_us(1:end-1 , idx_trial) * beta_drift;
sim_params_vec = [dt, 1, 0, true_bound, 0];

rt_simulated = nan(n_sims, 1);
rt_mex = nan(n_sims, 1);

for idx_sims = 1:n_sims

    rt_simulated(idx_sims) = drift_diff_new('mu_t', sm_signal(idx_trial , 1:end - 1) * beta_drift, 'theta', true_bound, ...
        'dt', fps, 'T', points.censoring, 'ndt', true_ndt);

    rt_mex(idx_sims) = sim_ddm_seeded(mu_tv, sim_params_vec, 1, idx_sims) + true_ndt;

end

rt_mex(isnan(rt_mex)) = 11;
rt_simulated(isnan(rt_simulated)) = 11;
rt_simulated(rt_simulated < points.truncation) = nan;
rt_mex(rt_mex < points.truncation) = nan;

bin_size = 1/60;
figure('Color', 'w', 'Position', [10 10 610 210])
hold on
histogram(rt_simulated(~isnan(rt_simulated)), time_vector(1:end-1) - dt/2, ...
    'Normalization', 'pdf', 'FaceColor', 'r', 'EdgeColor', 'none', ...
    'FaceAlpha', 0.3);
histogram(rt_mex(~isnan(rt_mex)), time_vector(1:end-1) - dt/2, ...
    'Normalization', 'pdf', 'FaceColor', 'b', 'EdgeColor', 'none', ...
    'FaceAlpha', 0.3);

params = param_res(true_bound, 'points', points, 'dt', dt);
[pde_results] = ddm_pdf_from_trace([0, true_bound, true_ndt], mu_tv, params);

plot(pde_results.t, pde_results.pdf, 'k--', 'LineWidth', 1);
stem(pde_results.t(end) + fps, pde_results.survival, 'k')

plot(time_vector_us(1:end-1), mu_tv, 'k--')
apply_generic(gca, 'no_y', true)

% Now we will repeat the analysis but with a model with stationary evidence

mu_st =  bouts_proc.sm(idx_trial) * ones(size(mu_tv));

for idx_sims = 1:n_sims

    rt_simulated(idx_sims) = drift_diff_new('mu_t', mu_st, 'theta', true_bound, ...
        'dt', dt, 'T', points.censoring, 'ndt', true_ndt);

    rt_mex(idx_sims) = sim_ddm_seeded(mu_st, sim_params_vec, 1, idx_sims) + true_ndt;

end

rt_mex(isnan(rt_mex)) = 11;
rt_simulated(isnan(rt_simulated)) = 11;
rt_simulated(rt_simulated < points.truncation) = nan;
rt_mex(rt_mex < points.truncation) = nan;

figure('Color', 'w', 'Position', [10 10 610 210])
hold on
histogram(rt_simulated(~isnan(rt_simulated)), time_vector_us(1:end-1) - dt/2, ...
    'Normalization', 'pdf', 'FaceColor', 'r', 'EdgeColor', 'none', ...
    'FaceAlpha', 0.3);
histogram(rt_mex(~isnan(rt_mex)), time_vector_us(1:end-1) - dt/2, ...
    'Normalization', 'pdf', 'FaceColor', 'b', 'EdgeColor', 'none', ...
    'FaceAlpha', 0.3);

params = param_res(true_bound, 'points', points, 'dt', dt);
[pde_results] = ddm_pdf_from_trace([0, true_bound, true_ndt], mu_st, params);

plot(pde_results.t, pde_results.pdf, 'k--', 'LineWidth', 1)
stem(pde_results.t(end) + fps, pde_results.survival, 'k')
trapz(pde_results.t(pde_results.t >= points.truncation & pde_results.t <= points.censoring), pde_results.pdf(pde_results.t >= points.truncation & pde_results.t <= points.censoring)) +  pde_results.survival
pde_mass = sum(pde_results.pdf(pde_results.idx_trunc:end)) * dt + pde_results.survival

apply_generic(gca, 'no_y', true)

[~, f, fd] = nll_fly_ddm_newer([1, true_bound, true_ndt], bouts_proc(idx_trial, :), points, 'model_sddm6', 'iid', 'p', []);
trapz(fd(fd >= points.truncation & fd <= points.censoring), f(fd >= points.truncation & fd <= points.censoring)) +  f(end)

plot(fd, f,  'r.-', 'LineWidth', 1)

% Ensure same time alignment
f_pde = pde_results.pdf(pde_results.idx_trunc:end);
% Interpolate f_theory (from nll) onto the same grid if necessary
f_theory_interp = interp1(fd, f, pde_results.t(pde_results.idx_trunc:end));

% Calculate difference
error_norm = sum((f_pde - f_theory_interp).^2);
fprintf('Difference between methods: %e\n', error_norm);
