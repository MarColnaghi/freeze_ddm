function [lls] = evaluate_likelihood_parallel(synth_data, motion_cache, ddm_params, model_results, select_indexes)

col = cmapper;
headers = {'id', 'll_tv', 'll_st', 'll_cf'};
model_str = sprintf('model_%s', model_results.fitted_model);
model = eval(model_str);
estimates_mean = model_results.estimates_mean;
est_params = table2array(estimates_mean(:, find(~ismissing(estimates_mean))));
ncomp_vars = evaluate_model(model, est_params, synth_data);
chunk_len = length(ddm_params.time_vector);
points = model_results.points;
points.truncation = [];

if isempty(select_indexes)
    trial_indexes = 1:ddm_params.eval_trials;
else
    trial_indexes = select_indexes;
end
n_trials = length(trial_indexes);

ll_id = cell(n_trials, 1);
ll_tv = zeros(n_trials, 1);
ll_st = zeros(n_trials, 1);
ll_cf = zeros(n_trials, 1);

% Progress tracking setup
progress = 0;
dq = parallel.pool.DataQueue;
afterEach(dq, @updateProgress);

% Run in parallel
parfor i = 1:n_trials
    idx_trials = trial_indexes(i);
    trial_data = synth_data(idx_trials, :);
    ons = trial_data.onsets;
    motion_ts = motion_cache(trial_data.fly);
    sm_raw = motion_ts(ons:ons + chunk_len - 1);
    tndt_s = ncomp_vars.tndt(idx_trials);
    rts_tv = nan(ddm_params.num_sims, 1);
    rts_st = nan(ddm_params.num_sims, 1);

    for j = 1:ddm_params.num_sims
        if rand < ncomp_vars.pmix(idx_trials)
            mu_t = estimates_mean.mu1_sm .* sm_raw;
            mu_st = estimates_mean.mu1_sm .* trial_data.sm;
            theta_s = ncomp_vars.theta1(idx_trials);
        else
            mu_t = estimates_mean.mu2_sm .* sm_raw;
            mu_st = estimates_mean.mu2_sm .* trial_data.sm;
            theta_s = ncomp_vars.theta2(idx_trials);
        end
        rts_st(j) = drift_diff_new('mu', mu_st, 'theta', theta_s, 'z', ddm_params.z, ...
            'dt', ddm_params.dt, 'T', ddm_params.T, 'ndt', tndt_s);
        rts_tv(j) = drift_diff_new('mu_t', mu_t, 'theta', theta_s, 'z', ddm_params.z, ...
            'dt', ddm_params.dt, 'T', ddm_params.T, 'ndt', tndt_s);
    end

    ll_st_i = kde_logL_censored(rts_st, trial_data.durations_s, ddm_params.kde_grid);
    ll_tv_i = kde_logL_censored(rts_tv, trial_data.durations_s, ddm_params.kde_grid);
    [ll_cf_i, ~, ~] = nll_fly_ddm_newer(est_params, trial_data, points, model_str, 'iid', 'p', []);
    ll_cf_i = -ll_cf_i;

    % Store results
    ll_id{i} = trial_data.id;
    ll_tv(i) = ll_tv_i;
    ll_st(i) = ll_st_i;
    ll_cf(i) = ll_cf_i;

    send(dq, i); % send index for possible debugging, but not used here
end

% Build result table
lls = table(vertcat(ll_id{:}), ll_tv, ll_st, ll_cf, 'VariableNames', headers);

% Save final output
save('lls_checkpoint.mat', 'lls');
save('lls_final.mat', 'lls');
disp('Final results saved.');

% Nested progress update function
    function updateProgress(~)
        progress = progress + 1;
        if mod(progress, 5) == 0 || progress == n_trials
            fprintf('Progress: %d / %d\n', progress, n_trials);
            % Safe bounds in case progress > n_trials due to async updates
            upto = min(progress, n_trials);
            fprintf('Partial TV-ST diff: %.3f\n', sum(ll_tv(1:upto) - ll_st(1:upto), 'omitnan'));
            fprintf('TV better in %.1f%% trials\n', 100 * mean((ll_tv(1:upto) - ll_st(1:upto)) > 0));
        end
        if mod(progress, 200) == 0
            save('lls_checkpoint.mat', 'll_id', 'll_tv', 'll_st', 'll_cf');
            disp('Checkpoint saved.');
        end
    end
end
