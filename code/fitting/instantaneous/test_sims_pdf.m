
true_drift_scale = 1.5;
true_lambda = 0;
true_bound = 3;

% -- Likelihood / Signal Parameters (Coarse) --
fixed.dt = 1/60;       % 1ms for signal and PDE
fixed.dx = 0.01;
fixed.sigma_sq = 1.0;
%     fixed.bound = 1.0;
fixed.x0 = 0.0;

% -- Simulation Parameters (Fine) --
sim.dt = 0.00005;       % 0.05ms for ground truth simulation

% PDE Grid Setup
x_min = -7;
grid_size = round((true_bound - x_min) / fixed.dx) + 1;
start_idx = round((fixed.x0 - x_min) / fixed.dx) + 1;

% Save for objective function
fixed.x_min = x_min;
fixed.grid_size = grid_size;
fixed.start_idx = start_idx;

% Truncation and T_max
fixed.T_trunc = 0;
fixed.T_max = 10.5 - 1/60;
bouts_le = bouts_le(bouts_le.durations_s > 0.5, :);

%%
paths = path_generator('folder', 'fitting_freezes/le');
motion_cache = importdata(fullfile(paths.cache_path, 'motion_cache.mat'));


chunk_len = 630;
for idx_trials = 1:height(bouts_le)

    trial_data = bouts_le(idx_trials, :);
    ons = trial_data.onsets;
    off = trial_data.ends;

    motion_ts = motion_cache(trial_data.fly) ./ 10;
    sm_raw{idx_trials} = motion_ts(ons:ons + chunk_len - 1);

end

%%
p = [1.1, 00, 2.6];
i = 2;
out = ddm_pdf_from_trace(p, sm_raw{i}, fixed);

figure('color','w'); hold on; plot(out.t, out.pdf); xlabel('t'); ylabel('p(T=t)'); plot(out.t, sm_raw{i});

title('First-passage time PDF (renormalized)');

mu_tv = p(1) *  sm_raw{i};
 
rt = nan(10000, 1);
for idx_sims = 1:1000000
    rt(idx_sims) = drift_diff_new('mu_t', mu_tv(1:end-1), 'theta', p(3), 'z', 0, ...
            'dt', 1/60, 'T', max(out.t), 'ndt', 0);
end

histogram(rt, -1/120:3/60:11, 'Normalization', 'pdf')

rt(isnan(rt)) = 13;
[kde, xkde] = ksdensity(rt, 'BandWidth', 0.15); %, 'censoring', isnan(rt));

trapz(out.t, out.pdf)
trapz(xkde(xkde < 11), kde(xkde < 11))
plot(xkde, kde)

%%
figure('color','w'); 
for idx = 20 %:5:50
    p = [4.2, 0, .3];
    i = 33;
    out = ddm_pdf_from_trace(p, sm_raw{i}, fixed);

    hold on; plot(out.t, out.pdf); xlabel('t'); ylabel('p(T=t)'); plot(out.t, sm_raw{i});
    trapz(out.t(~isnan(out.pdf)), out.pdf(~isnan(out.pdf)));

    title('First-passage time PDF (renormalized)');
end
mu_tv = p(1) *  sm_raw{i};

rt = nan(10000, 1);
for idx_sims = 1:50000
    rt(idx_sims) = extrema_detection_new('mu_t', mu_tv(1:end-1), 'theta', p(3), 'z', 0, ...
            'dt', 1/60, 'T', max(out.t), 'ndt', 0);
end


histogram(rt, -1/120:1/60:11, 'Normalization', 'pdf')
ec.soc_mot_array = sm_raw{i}';
points.censoring = max(out.t);
points.truncation = 0;
bouts_le.sm = bouts_le.avg_sm_freeze_norm;
bouts_le.fs = bouts_le.avg_fs_1s_norm;
bouts_le.ln = bouts_le.nloom_norm;
bouts_le.ls = bouts_le.sloom_norm;
bouts_le.intercept = ones(height(bouts_le),1);
bouts_le.smp = bouts_le.avg_sm_freeze_norm;
bouts_le.onsets = bouts_le.onsets;
bouts_le.fly = bouts_le.fly;
[~, f, fd] =  nll_fly_ddm_newer([p(1) p(3) 0], bouts_le(i, :), points, 'model_ed1', 'iid', 'p', ec);

plot(fd, f * 60)

%%
figure('color','w'); 
for idx = 20 %:5:50
    p = [0.9, 0, 1.6];
    i = 33;
    sm_raw{i} = ones(1, length(sm_raw{i}));
    out = ddm_pdf_from_trace(p, sm_raw{i}, fixed);

    hold on; plot(out.t, out.pdf); xlabel('t'); ylabel('p(T=t)'); plot(out.t, sm_raw{i});
    trapz(out.t(~isnan(out.pdf)), out.pdf(~isnan(out.pdf)));

    title('First-passage time PDF (renormalized)');
end
mu_tv = p(1) *  sm_raw{i};

for idx_sims = 1:100000
    rt(idx_sims) = drift_diff_new('mu_t', mu_tv(1:end-1), 'theta', p(3), 'z', 0, ...
            'dt', 1/60, 'T', max(out.t), 'ndt', 0);
end

histogram(rt, -1/120:3/60:11, 'Normalization', 'pdf')

[~, f, fd] =  nll_fly_ddm_newer([p(1) p(3) 0], bouts_le(i, :), points, 'model_sddm0', 'iid', 'p', ec);

plot(fd + 1/120, f)

%%
rt = nan(10000, 1);
for idx_sims = 1:50000
    rt(idx_sims) = extrema_detection_new('mu_t', mu_tv(1:end-1), 'theta', p(3), 'z', 0, ...
            'dt', 1/60, 'T', max(out.t), 'ndt', 0);
end


histogram(rt, -1/120:1/60:11, 'Normalization', 'pdf')
ec.soc_mot_array = sm_raw{i}';
points.censoring = max(out.t);
points.truncation = 0;
bouts_le.sm = bouts_le.avg_sm_freeze_norm;
bouts_le.fs = bouts_le.avg_fs_1s_norm;
bouts_le.ln = bouts_le.nloom_norm;
bouts_le.ls = bouts_le.sloom_norm;
bouts_le.intercept = ones(height(bouts_le),1);
bouts_le.smp = bouts_le.avg_sm_freeze_norm;
bouts_le.onsets = bouts_le.onsets;
bouts_le.fly = bouts_le.fly;
[~, f, fd] =  nll_fly_ddm_newer([p(1) p(3) 0], bouts_le(i, :), points, 'model_ed1', 'iid', 'p', ec);

plot(fd, f * 60)
