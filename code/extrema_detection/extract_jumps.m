clear all
close all
idx_seed = randi(100);

exporting = true;
col = cmapper();
sim_params.rng = idx_seed;
rng(sim_params.rng);

% Model
model = 'dddm2';
select_run = 'run03';
gen_data = 'fr';

% Load important files
paths = path_generator('folder', fullfile('/fitting_freezes/le', model, select_run));
load(fullfile(paths.results, 'surrogate.mat'))
load('/Users/marcocolnaghi/PhD/freeze_ddm/model_results/momentary_integration/surrogate/dddm2/run13/sims_tv/y.mat')

% Load the motion ts
sim_params.motion_cache_path = fullfile(paths.dataset, 'motion_cache.mat');
load(sim_params.motion_cache_path)

d = 240;
n_selected_comparisons = 50;
y = y(y.durations_s > d/60 & y.durations_s <= 10.5, :);
sm_cell = extract_sm_cellarray(y, motion_cache);
col = cbrewer2('Spectral', n_selected_comparisons);

for idx_bout = 10:40

    sm = motion_cache(y.fly(idx_bout));
    sm_bout = sm(y.onsets(idx_bout):y.onsets(idx_bout) + round(y.durations_s(idx_bout) .* 60));
    offset = length(sm_bout); onset = offset - d;
    v1 = sm_bout(onset:offset);
    similarity_value = nan(1, height(y));
    similarity_sort = nan(1, height(y));

    for idx_comparison = 1:height(y)
        extracted_sm = sm_cell{idx_comparison};

        if length(extracted_sm) < d
            similarity_value(idx_comparison) = nan;
            similarity_sort(idx_comparison) = nan;
            disp('hello')
        else

            similarity_framebframe = nan(1, length(extracted_sm) - d);

            for idx_frame = 1:length(extracted_sm) - d
                v2 = extracted_sm(idx_frame:idx_frame + d);
                similarity_framebframe(idx_frame) = norm(v1 - v2);
            end

            [closest_similarity, best_frame] = min(similarity_framebframe);
            similarity_value(idx_comparison) = closest_similarity;
            similarity_sort(idx_comparison) = best_frame;

        end

    end

    [sorted_similarity_magnitude, sorted_similarity_idx] = sort(similarity_value);
    best_similarities = sorted_similarity_idx(1:n_selected_comparisons);

    figure
    tiledlayout(2,1)
    nexttile
    hold on

    for idx_tops = 1:length(best_similarities)
        extracted_sm = sm_cell{best_similarities(idx_tops)};
        %v2 = extracted_sm(similarity_sort(best_similarities(idx_tops)):similarity_sort(best_similarities(idx_tops)) + d);
        %plot(1:length(v2), v2, 'Color', col(idx_tops, :))
        init = similarity_sort(best_similarities(idx_tops));
        v2 = extracted_sm(similarity_sort(best_similarities(idx_tops)):similarity_sort(best_similarities(idx_tops)) + d);
        plot(-init + 1: - init + length(extracted_sm), extracted_sm, 'Color', col(idx_tops, :))

        initial_bump(idx_bout, idx_tops) = sum(extracted_sm(1:similarity_sort(best_similarities(idx_tops))));
    end

    plot(1:length(v1), v1, 'k--', 'LineWidth', 2)
    nexttile
    scatter(initial_bump(idx_bout, :), y.durations_s(best_similarities,:), 60, col, 'filled')

end

%%
figure
hold on
histogram(similarity_value, 0:0.1:50)





%%
figure 
hold on
for n_mov = 0:4 
    y = surrogate(surrogate.moving_flies == n_mov, :);
    i = 0;
    for idx_bout = 1:height(y)
        i = i + 1;
        sm = motion_cache(y.fly(idx_bout));
        sm_bout(i, :) = sm(y.onsets(idx_bout):y.onsets(idx_bout) + y.durations(idx_bout));
    end
    plot(mean(sm_bout))

end

