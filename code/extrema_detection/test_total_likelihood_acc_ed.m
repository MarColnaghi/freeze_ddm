% In this script we have a look at individual freezes and compute their
% likelihood under the two alternative models

% Load the table first. We will take advantage of an already existing
% dataset.
threshold_imm = 2; threshold_mob = 2; threshold_pc = 4; id_code = sprintf('imm%d_mob%d_pc%d', threshold_imm, threshold_mob, threshold_pc);

col = cmapper();

color_list = {col.timevarying_sm, col.extremadetection};
model_list = {'ded2', 'dddm2' };
run_list = {'run10', 'run13' };
label_list = {'Extrema Detection','Integration'};

paths = path_generator('folder', 'fitting_freezes/le');

results_folder_1 = fullfile(paths.results, model_list{1}, run_list{1});
results_folder_2 = fullfile(paths.results, model_list{2}, run_list{2});

freezes_1 = importdata(fullfile(results_folder, 'surrogate.mat'));
freezes_2 = importdata(fullfile(results_folder, 'surrogate.mat'));

if isequaln(freezes_1, freezes_2)
    fprintf('The two freeze tables match \n')
else
    fprintf("The two freeze tables don't match \n")
    return

end


for idx_model = 1:length(model_list)

    results_folder = fullfile(paths.results, model_list{idx_model}, run_list{idx_model});

    extra = importdata(fullfile(results_folder, 'extra.mat'));

    model_results = importdata(fullfile(results_folder, sprintf('fit_results_%s.mat', model_list{idx_model})));
    freezes = importdata(fullfile(results_folder, 'surrogate.mat'));
    freezes.sm = extra.soc_mot_array;
    freezes.durations(freezes.durations_s >= model_results.points.censoring) = 631;

    est_params = table2array(model_results.estimates_mean(:, ~ismissing(model_results.estimates_mean)));

    model_func = str2func(strcat('model_', model_list{idx_model}));
    model = model_func();
    out = evaluate_model(model, model_results.estimates_mean, freezes);

    ll = nan(height(freezes), 1);
    nll_st = nan(height(freezes), 1);

    for idx_bout = 1:height(freezes)

        out_f = out(idx_bout, :);
        out_f.tndt = 0;
        fprintf('bout: %d \n', idx_bout)
        freeze = freezes(idx_bout, :);
            
        ec.soc_mot_array = extra.soc_mot_array(idx_bout,:);
        is_censored = freeze.durations_s >= model_results.points.censoring;

        if contains(model_list{idx_model}, 'ddm')

            [pdf] = compute_pdf_tv_ddm(out_f);
            [ll(idx_bout)] = compute_likelihood(pdf, is_censored, freeze.durations_s);

            model_results.points.truncation = 0;

            % Check if ndt = 0 works better: it does.

            % est_params(end) = 0;
            [~, f, fd] = nll_fly_ddm_newer(est_params, freeze, model_results.points, strcat('model_', model_list{idx_model}), 'iid', 'p', ec);
            [nll_st(idx_bout)] = nll_fly_ddm_newer(est_params, freeze, model_results.points, strcat('model_', model_list{idx_model}), 'iid', '', ec);

            fprintf('ll difference: %d \n', sum(ll + nll_st, 1, 'omitnan'))
 
%             figure
%             plot(fd, f)
%             hold on
%             plot(pdf.t, pdf.ddm)
%             xline(freezes(idx_bout, :).durations_s)
%             ll(idx_bout) + nll_st(idx_bout)
%             hold off

%             plot(pdf.t, pdf.ddm)
%             hold on
%             plot(fd, f)

            % Check that the pdf matches the theoretical one: it does! be
            % aware of the logit function transformation.

%             mu1 = 1.4;
%             mu2 = 2;
%             out_f.mu1 = mu1 * ones(1,630);
%             out_f.mu2 = mu2 * ones(1,630);
% 
%             [pdf] = compute_pdf_tv_ddm(out_f);
%             [ll(idx_bout)] = compute_likelihood(pdf, is_censored, freeze.durations_s);
%             [~, f, fd] = nll_fly_ddm_newer([mu1, out_f.theta1, mu2, out_f.theta2, log(out_f.pmix / (1 - out_f.pmix)), out_f.tndt], freezes(idx_bout, :), model_results.points, 'model_dddm0', 'iid', 'p', ec);
% 
%             figure
%             hold on
%             plot(pdf.t, pdf.ddm)
%             plot(fd, f)
% 
%             n_sims = 200000;
%             coin_toss = rand(n_sims, 1);
%             rt = nan(n_sims, 1);
%             for idx_sims = 1:n_sims
% 
%                 if coin_toss(idx_sims) < out_f.pmix
%                     rt(idx_sims) = drift_diff_new('mu_t', out_f.mu1, 'theta', out_f.theta1, 'z', 0, ...
%                         'dt', 1/60, 'T', 10.5, 'ndt', out_f.tndt);
%                 else
%                     rt(idx_sims) = drift_diff_new('mu_t', out_f.mu2, 'theta', out_f.theta2, 'z', 0, ...
%                         'dt', 1/60, 'T', 10.5, 'ndt', out_f.tndt);
%                 end
% 
%             end
% 
%             histogram(rt, -1/120:1/60:11, 'Normalization', 'pdf')

        elseif contains(model_list{idx_model}, 'ed')

            model_results.points.truncation = [];

            [~, f, fd] = nll_fly_ddm_newer(est_params, freezes(idx_bout, :), model_results.points, 'model_ded2', 'iid', 'p', ec);
            [nll_st(idx_bout)] = nll_fly_ddm_newer(est_params, freeze, model_results.points, strcat('model_', model_list{idx_model}), 'iid', '', ec);
            trapz(fd(1:end-1),f(1:end-1)) + f(end)
            figure
            plot(fd, f)
            hold on
            %plot(fd, ec.soc_mot_array)

        end

    end
end


