clear all
close all
% Load the table first. We will take advantage of an already existing
% dataset.

col = cmapper();

threshold_imm = 3; threshold_mob = 3; threshold_pc = 4; id_code = sprintf('imm%d_mob%d_pc%d', threshold_imm, threshold_mob, threshold_pc);
paths = path_generator('folder', 'fitting_tests/sddm', 'bouts_id', id_code);
load(fullfile(paths.dataset, 'bouts.mat'));
bouts_proc = data_parser_new(bouts, 'type', 'immobility', 'period', 'loom', 'window', 'le');
bouts_proc = bouts_proc(bouts_proc.nloom < 19, :);

% Load the motion ts
sim_params.motion_cache_path = fullfile(paths.dataset, 'motion_cache.mat');
load(sim_params.motion_cache_path)

% Only save useful variables in the table
y = table;
y.sm = bouts_proc.avg_sm_freeze_norm;
y.smp = bouts_proc.avg_sm_freeze_norm;
y.fs = bouts_proc.avg_fs_1s_norm;
y.ln = bouts_proc.nloom_norm;
y.ls = bouts_proc.sloom_norm;
y.intercept = ones(height(y),1);   % N‑by‑1 column of ones
predictors = y.Properties.VariableNames;

link_linear = @(x) x;     % log link for bound height
link_logistic = @(x) 1./(1 + exp(-x));     % log link for bound height

% For mu
model.mu = struct( ...
    'predictors', {{ ...
    struct('name', 'sm'), ...
    struct('name', 'intercept') ...
    }}, ...
    'mask', [1 0 0 0 0 1], ...
    'ground_truth', [2 0.4], ...
    'link', link_linear ...
    );

% For theta
model.theta = struct( ...
    'predictors', {{ ...
    struct('name', 'fs'), ...
    struct('name', 'ls'), ...
    struct('name', 'intercept') ...
    }}, ...
    'mask', [0 0 1 0 1 1], ...
    'ground_truth', [1 1 2], ...
    'link', link_linear ...
    );

% Non decision time
model.tndt = struct( ...
    'predictors', {{ struct('name', 'intercept') }}, ...
    'mask', [0 0 0 0 0 1], ...
    'ground_truth', 0, ...
    'link', link_linear ...
    );

[gt, lbl] = get_ground_truth_vector(model);
x = gt(~isnan(gt));
ncomp_vars = evaluate_model(model, x, y);
gt_table = array2table(gt, 'VariableNames', lbl);

% Specify the seed
rng(1);

% General simulation parameters
sim_params.n_trials = height(bouts_proc);
sim_params.dt = 1/60;
sim_params.T = 30;
sim_params.time_vector = 0:sim_params.dt:sim_params.T;
sim_params.z = 0;
sim_params.theta_factor = 12.1;
sim_params.gt_table = gt_table;

% Simulation settings
sim_params.kde_grid = 0:1/600:120;
sim_params.eval_trials = sim_params.n_trials;
sim_params.num_sims = 20000;

% Censoring/Truncation
sim_params.points.truncation = [];
sim_params.points.censoring = sim_params.T;

% Initialize outputs
rt = table;
rt.ac = nan(sim_params.n_trials, 1);
rt.ed = nan(sim_params.n_trials, 1);
sm_raw = cell(sim_params.n_trials, 1);

% Extract Social Motion TimeSeries
chunk_len = length(sim_params.time_vector);

tic
for idx_trials = 1:height(bouts_proc)
%     figure
%     hold on
    ons = bouts_proc.onsets(idx_trials);

    sum_motion = motion_cache(bouts_proc.fly(idx_trials));
    sm_raw{idx_trials} = sum_motion(ons:ons + chunk_len - 1);
    sm_chunk = sm_raw{idx_trials};

    mu_tv = gt_table.mu_sm .* sm_chunk + gt_table.mu_intercept;
    theta_s = ncomp_vars.theta(idx_trials);
    tndt_s = ncomp_vars.tndt(idx_trials);

    % Simulate RT from full DDM
    [rt.ac(idx_trials), traj_st] = drift_diff_new('mu_t', mu_tv, 'theta', theta_s, ...
        'z', sim_params.z, 'dt', sim_params.dt, 'T', sim_params.T, 'ndt', tndt_s, 'seed', idx_trials, 'sigma', 1);

    [rt.ed(idx_trials), traj_ed] = extrema_detection_new('mu_t', mu_tv , 'theta', theta_s / sim_params.theta_factor, ...
        'z', sim_params.z, 'dt', sim_params.dt, 'T', sim_params.T, 'ndt', tndt_s, 'seed', idx_trials, 'sigma', 1);

end
toc

rt = [rt ncomp_vars];
fh = figure('color', 'w', 'Position', [100, 100, 600, 600]);

hold on
histogram(rt.ed, sim_params.dt/2:sim_params.dt * 5 :sim_params.T + 1, 'FaceColor', col.extremadetection, 'EdgeColor', 'none')
histogram(rt.ac, sim_params.dt/2:sim_params.dt * 5 :sim_params.T + 1, 'FaceColor', col.timevarying_sm, 'EdgeColor', 'none')
xline(mean(rt.ed, 'omitnan'), 'Color', col.timevarying_sm, 'LineWidth', 2)
xline(mean(rt.ac, 'omitnan'), 'Color', col.extremadetection, 'LineWidth', 2)
xlabel('Duration')
ylabel('Count')
apply_generic(gca, 24);

%%
rt.ac(isnan(rt.ac)) = sim_params.T + 1;
d = 180;
frames_2b_exp = 60;
n_selected_comparison = 120;
extra_frames = d - frames_2b_exp - 1;

rt_temp = rt(rt.ac / sim_params.dt > d & rt.ac / sim_params.dt < sim_params.T * 60, :);
freezes_durations = rt_temp.ac / sim_params.dt;

for idx_bout = 1:height(rt_temp)

    sm = sm_raw{idx_bout};
    freeze_duration = freezes_durations(idx_bout);
    sm_bout = sm(1:freeze_duration);

    template_onset = length(sm_bout) - d + 1;

    v1 = sm_bout(template_onset:end);

    similarity_value = nan(1, height(rt_temp));
    similarity_sort = nan(1, height(rt_temp));

    all_cropped = cell(1, height(rt_temp));

    for idx_comparison = 1:height(rt_temp)

        comparison_sm = sm_raw{idx_comparison};

        start = 1;
        ending = freezes_durations(idx_comparison) + extra_frames;
        comparison_sm_cropped = comparison_sm(start:ending);

        dots = conv(comparison_sm_cropped, flipud(v1), 'valid');
        ssq  = conv(comparison_sm_cropped.^2, ones(d, 1), 'valid');
        v1sq = sum(v1.^2);

        similarity_framebframe_vectorized = sqrt(max(v1sq + ssq - 2*dots, 0));

        [closest_similarity, best_frame] = min(similarity_framebframe_vectorized);
        similarity_value(idx_comparison) = closest_similarity;
        similarity_sort(idx_comparison) = best_frame;
        all_cropped{idx_comparison} = comparison_sm_cropped;

    end

    [~, best_similarities] = mink(similarity_value, n_selected_comparisons);


    for idx_tops = 1:length(best_similarities)

        comparison_sm_cropped = all_cropped{best_similarities(idx_tops)};

        starting_frame = similarity_sort(best_similarities(idx_tops));

        initial_bump(idx_bout, idx_tops) = sum(comparison_sm_cropped(1:starting_frame - 1));

        ending_duration(idx_bout, idx_tops) = length(comparison_sm_cropped(1 : freezes_durations(best_similarities(idx_tops)))) - length(comparison_sm_cropped(1:starting_frame - 1));

    end

    [init_bump_sorted, init_bump_idx] = sort(initial_bump(idx_bout, :), 'descend');

    corr_bout = corrcoef(init_bump_sorted, ending_duration(idx_bout, init_bump_idx) ./60);

    ending_duration_sorted(idx_bout, :) =  ending_duration(idx_bout, init_bump_idx);

    corr_array_surr(idx_bout) = corr_bout(1,2);
end

