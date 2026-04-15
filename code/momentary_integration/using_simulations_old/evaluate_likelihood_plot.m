function [lls, ax] = evaluate_likelihood_plot(synth_data, motion_cache, ddm_params, model_results, select_indexes, varargin)

opt = inputParser;
addParameter(opt, 'export', false);
addParameter(opt, 'legend', true);

parse(opt, varargin{:});

export = opt.Results.export;
legend = opt.Results.legend;


col = cmapper;
lls = table;
headers = {'id','ll_tv', 'll_st', 'll_cf'};


model_str = sprintf('model_%s', model_results.fitted_model);
model = eval(model_str);
estimates_mean = model_results.estimates_mean;
est_params = table2array(estimates_mean(:, find(~ismissing(estimates_mean))));
ncomp_vars = evaluate_model(model, est_params, synth_data);
chunk_len = length(ddm_params.time_vector) - 1;

points = model_results.points;
points.truncation = [];

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
    %histogram(rts_st,  1/120:2/60:points.censoring, 'DisplayStyle', 'stairs', 'Normalization', 'pdf', 'LineWidth', 1, 'EdgeColor', col.stationary_sm)
    %histogram(rts_tv,  1/120:2/60:points.censoring, 'DisplayStyle', 'stairs', 'Normalization', 'pdf', 'LineWidth', 1, 'EdgeColor', col.timevarying_sm)
    
    %plot(fd, f, 'k--')
    xlabel('Duration (s)')
    ylabel('density')

    ax = gca;
    apply_generic(ax)

    create_vector_for_imagesc = nan(size(ddm_params.kde_grid));
    create_vector_for_imagesc(1:length(sm_raw)) = sm_raw;
    imagesc(ddm_params.kde_grid, -0.55, create_vector_for_imagesc, [0 2]);  % set XData = time_vector
    colormap(col.vars.sm)

    xlim([0 11])
    ylim([-0.3 2.0])
    xticks([0 10])

    ax.Clipping= 'on';
    set(ax ,'Layer', 'Top')
    ax.YAxis.Visible = 'off';
    ax.XLabel.Position(2) = -0.45;

    if ll_tv - ll_st > 0
        text(trial_data.durations_s + 0.2, 1.5, sprintf('$LL_{st}: %.2f$', ll_st), 'HorizontalAlignment', 'left', 'FontSize', 20, 'Interpreter', 'latex', 'Color',  'k')
        text(trial_data.durations_s + 0.2, 1.75, sprintf('$LL_{tv}: %.2f$', ll_tv), 'HorizontalAlignment', 'left', 'FontSize', 20, 'Interpreter', 'latex', 'Color',  col.timevarying_sm)

    else
        text(trial_data.durations_s+ 0.2, 1.5, sprintf('$LL_{st}: %.2f$', ll_st), 'HorizontalAlignment', 'left', 'FontSize', 20, 'Interpreter', 'latex', 'Color',  col.stationary_sm)
        text(trial_data.durations_s+ 0.2, 1.75, sprintf('$LL_{tv}: %.2f$', ll_tv), 'HorizontalAlignment', 'left', 'FontSize', 20, 'Interpreter', 'latex', 'Color',  'k')
    end

    xline(trial_data.durations_s, 'LineWidth', 2, 'LineStyle', '--', 'Color', 'k')
    scatter(trial_data.durations_s, ax.YLim(1) - 0.15, 100, 'filled', '^', 'Clipping', 'off', 'MarkerFaceColor', 'k')
    text(trial_data.durations_s + 0.5, ax.YLim(1) - 0.15, sprintf('id: %d', trial_data.id))

    if legend
        xticklabels([0 10])
    else
        xticklabels({''})
        xlabel('')
    end


    drawnow



end
end
