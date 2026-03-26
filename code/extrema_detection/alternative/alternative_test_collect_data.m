function data = alternative_test_collect_data(model_list, run_list)

paths_fit      = path_generator('folder', 'fitting_freezes/le_new');

% --- 2. Data Pre-loading & Preparation ---
data = struct();



for m = 1:length(model_list)
    
    data(m).dt = 1/60;
    data(m).T = 10.5;
    data(m).t = 0:data(m).dt:data(m).T;

    fps = 1/60;
    time_vector = 0:fps:data(m).T;
    n_steps = numel(data(m).t);

    res_folder = fullfile(paths_fit.results, model_list{m}, run_list{m});

    extra   = importdata(fullfile(res_folder, 'extra.mat'));
    res     = importdata(fullfile(res_folder, 'model_results.mat'));
    freezes = importdata(fullfile(res_folder, 'freeze.mat'));

    % --- STEP A: Check equality BEFORE modification ---
    data(m).raw_table = freezes;

    % --- STEP B: Now apply model-specific modifications ---
    freezes.sm_stat = freezes.sm;
    freezes.sm = freezes.sm * ones(1, n_steps);

    data(m).results = res;
    data(m).table   = freezes;
    data(m).func    = str2func(strcat('model_', model_list{m}));
    data(m).est     = table2array(res.estimates_mean(:, ~ismissing(res.estimates_mean)));
    data(m).extra_st   = freezes.sm;

    data(m).out_all_st = evaluate_model(data(m).func(), res.estimates_mean, freezes);

    sm_signal_us = interp1(time_vector(:), extra.soc_mot_array', data(m).t(:), 'nearest');

    % Use the specific censoring point for this specific model run
    freezes.durations_s(freezes.durations_s > res.points.censoring) = 660;

    data(m).extra_tv   = sm_signal_us';
    data(m).extra_tv_n = extra.soc_mot_array;
    freezes.sm = sm_signal_us';
    data(m).out_all_tv = evaluate_model(data(m).func(), res.estimates_mean, freezes);

end

end