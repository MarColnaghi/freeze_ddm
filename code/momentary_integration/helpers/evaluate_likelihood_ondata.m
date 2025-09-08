function [lls, ax] = evaluate_likelihood_ondata(synth_data, motion_cache, ddm_params, model_results, select_indexes)
    
col = cmapper;
lls = table;  % Make sure this is initialized
headers = {'id','ll_tv', 'll_st', 'll_cf'};

% progress = 0;
% dq = parallel.pool.DataQueue;
% afterEach(dq, @(~) updateProgress());
% 
% function updateProgress()
%     progress = progress + 1;
%     fprintf('Progress: %d / %d\n', progress, n_trials);
% end

model_str = sprintf('model_%s', model_results.fitted_model);
model = eval(model_str);
estimates_mean = model_results.estimates_mean;
est_params = table2array(estimates_mean(:, find(~ismissing(estimates_mean))));
ncomp_vars = evaluate_model(model, est_params, synth_data);
chunk_len = length(ddm_params.time_vector);

points = model_results.points;
points.truncation = [];

if isempty(select_indexes)

    for idx_trials = 1:ddm_params.eval_trials

        trial_data = synth_data(idx_trials, :);

        ons = trial_data.onsets;
        motion_ts = motion_cache(trial_data.fly);
        sm_raw = motion_ts(ons:ons + chunk_len - 1);

        % Determine model 1 or 2 based on pmix

        tndt_s = ncomp_vars.tndt(idx_trials);

        % Simulate RTs
        rts_tv = nan(ddm_params.num_sims,1);
        rts_st = nan(ddm_params.num_sims,1);
        
        for j = 1:ddm_params.num_sims

            if rand < ncomp_vars.pmix(idx_trials)
                mu_t = estimates_mean.mu1_sm .* sm_raw;
                mu_st = estimates_mean.mu1_sm .* trial_data.sm;
                theta_s = ncomp_vars.theta1(idx_trials);

            else
                mu_t = estimates_mean.mu2_sm  .* sm_raw;
                mu_st = estimates_mean.mu2_sm .* trial_data.sm;
                theta_s = ncomp_vars.theta2(idx_trials);

            end

            rts_st(j) = drift_diff_new('mu', mu_st, 'theta', theta_s, 'z', ddm_params.z, ...
                'dt', ddm_params.dt, 'T', ddm_params.T, 'ndt', tndt_s);
            rts_tv(j) = drift_diff_new('mu_t', mu_t, 'theta', theta_s, 'z', ddm_params.z, ...
                'dt', ddm_params.dt, 'T', ddm_params.T, 'ndt', tndt_s);
        end

        % Compute log-likelihoods
        [ll_tv, ~] = kde_logL_censored(rts_tv, trial_data.durations_s, ddm_params.kde_grid);
        [ll_st, ~] = kde_logL_censored(rts_st, trial_data.durations_s, ddm_params.kde_grid);
        [ll_cf, ~, ~] = nll_fly_ddm_newer(est_params, trial_data, points, model_str, 'iid', 'p', []);
        ll_cf = -ll_cf;

        lls(idx_trials, :) = table(trial_data.id, ll_tv, ll_st, ll_cf, 'VariableNames', headers);

        if mod(idx_trials, 5) == 0
            fprintf('freeze #%d \n', idx_trials)
            fprintf('Difference between tv and st is: %.3f \n', sum(lls{:, 2}, 1, 'omitnan') - sum(lls{:, 3}, 1, 'omitnan'))
            fprintf('Percentage in favour of tv: %.3f \n', sum(lls{:, 2} - lls{:, 3} > 0) ./ idx_trials);

        end
        if mod(idx_trials, 200) == 0
            save('lls.mat', 'lls')
            disp('saving... \n')

        end


    end

else
    for idx_trials = select_indexes'

        trial_data = synth_data(idx_trials, :);

        ons = trial_data.onsets;
        motion_ts = motion_cache(trial_data.fly);
        sm_raw = motion_ts(ons:ons + chunk_len - 1);

        % Determine model 1 or 2 based on pmix

        tndt_s = ncomp_vars.tndt(idx_trials);

        % Simulate RTs
        rts_tv = nan(ddm_params.num_sims,1);
        rts_st = nan(ddm_params.num_sims,1);
        
        for j = 1:ddm_params.num_sims

            if rand < ncomp_vars.pmix(idx_trials)
                mu_t = estimates_mean.mu1_sm .* sm_raw;
                mu_st = estimates_mean.mu1_sm .* trial_data.sm;
                theta_s = ncomp_vars.theta1(idx_trials);

            else
                mu_t = estimates_mean.mu2_sm  .* sm_raw;
                mu_st = estimates_mean.mu2_sm .* trial_data.sm;
                theta_s = ncomp_vars.theta2(idx_trials);

            end

            rts_st(j) = drift_diff_new('mu', mu_st, 'theta', theta_s, 'z', ddm_params.z, ...
                'dt', ddm_params.dt, 'T', ddm_params.T, 'ndt', tndt_s);
            rts_tv(j) = drift_diff_new('mu_t', mu_t, 'theta', theta_s, 'z', ddm_params.z, ...
                'dt', ddm_params.dt, 'T', ddm_params.T, 'ndt', tndt_s);
        end

        % Compute log-likelihoods
        [ll_tv, kde_tv] = kde_logL_censored(rts_tv, trial_data.durations_s, ddm_params.kde_grid);
        [ll_st, kde_st] = kde_logL_censored(rts_st, trial_data.durations_s, ddm_params.kde_grid);
        [ll_cf, f, fd] = nll_fly_ddm_newer(est_params, trial_data, points, model_str, 'iid', 'p', []);
        ll_cf = -ll_cf;

        lls(idx_trials, :) = table(trial_data.id, ll_tv, ll_st, ll_cf, 'VariableNames', headers);

      
        hold on

        plot(ddm_params.kde_grid, kde_st, 'LineWidth', 2, 'Color', col.stationary_sm)
        plot(ddm_params.kde_grid, kde_tv, 'LineWidth', 2, 'Color', col.timevarying_sm)
        histogram(rts_st,  1/120:2/60:points.censoring, 'DisplayStyle', 'stairs', 'Normalization', 'pdf', 'LineWidth', 1, 'EdgeColor', col.stationary_sm)
        histogram(rts_tv,  1/120:2/60:points.censoring, 'DisplayStyle', 'stairs', 'Normalization', 'pdf', 'LineWidth', 1, 'EdgeColor', col.timevarying_sm)
        plot(fd, f, 'k--')
        xlabel('Duration (s)')
        ylabel('density')

        xline(trial_data.durations_s, 'LineWidth', 2, 'LineStyle', '-.', ...
            'Label', {sprintf('LL_{st}: %.2f', ll_st), sprintf('LL_{tv}: %.2f', ll_tv)}, ...
            'LabelOrientation', 'horizontal', 'FontSize', 20, 'Interpreter', 'tex')
        ax(2) = gca;
        apply_generic(ax(2))

        create_vector_for_imagesc = nan(size(ddm_params.kde_grid));
        create_vector_for_imagesc(1:length(sm_raw)) = sm_raw;
        imagesc(ddm_params.kde_grid, -0.55, create_vector_for_imagesc, [0 2]);  % set XData = time_vector
        colormap(col.vars.sm)
        ax(1) = gca;
        ax(1).Clipping= 'on';
        set(ax ,'Layer', 'Top')

        apply_generic(ax(1))
        %             ax(1).XAxis.Visible = 'off';
        ax(1).YAxis.Visible = 'off';

        if ll_tv - ll_st > 0
            xline(trial_data.durations_s, 'LineWidth', 2, 'LineStyle', '--', 'Color', col.timevarying_sm)
        else
            xline(trial_data.durations_s, 'LineWidth', 2, 'LineStyle', '--', 'Color', col.stationary_sm)
        end

        linkaxes([ax(2), ax(1)], 'x')
        xlim([0 15])
        ylim([-0.1 1.6])

        drawnow
        xlim([0 15])
        ylim([-0.15 1.35])


        

    end

end

end