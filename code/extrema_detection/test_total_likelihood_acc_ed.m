% In this script we have a look at individual freezes and compute their
% likelihood under the two alternative models

% Load the table first. We will take advantage of an already existing
% dataset.
threshold_imm = 2; threshold_mob = 2; threshold_pc = 4; id_code = sprintf('imm%d_mob%d_pc%d', threshold_imm, threshold_mob, threshold_pc);

col = cmapper();

color_list = {col.timevarying_sm, col.extremadetection};
model_list = {'dddm2', 'ded2'};
run_list = {'run13', 'run10'};
label_list = {'Integration', 'Extrema Detection'};

paths = path_generator('folder', 'fitting_freezes/le');

results_folder = fullfile(paths.results, model_list{1}, run_list{1});
temp_results = importdata(fullfile(results_folder, sprintf('fit_results_%s.mat', model_list{1})));
surrogate = importdata(fullfile(results_folder, 'surrogate.mat'));


for idx_model = 1:length(model_list)

    results_folder = fullfile(paths.results, model_list{idx_model}, run_list{idx_model});

    model_results = importdata(fullfile(results_folder, sprintf('fit_results_%s.mat', model_list{idx_model})));
    freezes = importdata(fullfile(results_folder, 'surrogate.mat'));

    freezes.durations(freezes.durations_s >= model_results.points.censoring) = 631;

    extra = importdata(fullfile(results_folder, 'extra.mat'));

    est_params = table2array(model_results.estimates_mean(:, ~ismissing(model_results.estimates_mean)));
    ec = extra;

    model_func = str2func(strcat('model_', model_list{idx_model}));
    model = model_func();
    out = evaluate_model(model, model_results.estimates_mean, freezes);

    for idx_bout = 1:height(freezes)

        fprintf('bout: %d\n', idx_bout)
        freeze = freezes(idx_bout, :);
        ec.soc_mot_array = extra.soc_mot_array(idx_bout, :);
        soc_mot_vector = ec.soc_mot_array;

        is_censored = freeze.durations > model_results.points.censoring;

        if contains(model_list{idx_model}, 'ddm')

            compute_likelihood_tv_ddm(model_results, model, soc_mot_vector)

           

   

        end
        %
        %         l_ma(idx_bout) = nll_fly_ddm_newer([1 2 0], freeze, model_results.points, strcat('model_', 'ed0'), 'iid', '', ec);
        %         l_ju(idx_bout) = log_joint_likelihood_3D([1 1e-12 2], soc_mot_vector, freeze.durations_s, fixed);
        %
        %         [l_ma, l_ju]
    end
end


[sorted_ll, sorted_idx] = sort(diff(ll, [], 2));

%%

out = ddm_pdf_from_trace(p, sm_raw{i}, fixed);
