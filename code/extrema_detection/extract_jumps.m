clear all
close all

exporting = true;
col = cmapper();

% Model
model = 'dddm2';
select_run = 'run03';
gen_data = 'fr';

% Load important files
paths = path_generator('folder', fullfile('/fitting_freezes/le', model, select_run));
load(fullfile(paths.results, 'surrogate.mat'))
load('/Users/marcocolnaghi/PhD/freeze_ddm/model_results/momentary_integration/surrogate/dddm2/run11/sims_tv/y.mat')

% Load the motion ts
sim_params.motion_cache_path = fullfile(paths.dataset, 'motion_cache.mat');
load(sim_params.motion_cache_path)

% Parameters
d = 180; percentage_2b_exp = 1/3; extra_frames = round(percentage_2b_exp * d);
n_selected_comparisons = 120;
y = y(y.durations_s > d/60 & y.durations_s <= 10.5, :);
sm_cell = extract_sm_cellarray(y, motion_cache);
col = cbrewer2('Spectral', n_selected_comparisons);

% Preallocate once before the outer loop
initial_bump    = nan(height(y), n_selected_comparisons);
ending_duration = nan(height(y), n_selected_comparisons);
corr_array      = nan(height(y), 1);


for idx_bout = 1:height(y)

    sm = motion_cache(y.fly(idx_bout));
    sm_bout = sm(y.onsets(idx_bout):y.onsets(idx_bout) + round(y.durations_s(idx_bout) .* 60));
    
    offset = length(sm_bout); onset = offset - d;
    v1 = sm_bout(onset:offset);
    similarity_value = nan(1, height(y));
    similarity_sort = nan(1, height(y));
    all_cropped = cell(1, height(y));

    for idx_comparison = 1:height(y)
        comparison_sm = motion_cache(y.fly(idx_comparison));
        start = y.onsets(idx_comparison);
        ending = y.onsets(idx_comparison) + round(y.durations_s(idx_comparison) .* 60)  + extra_frames;
        comparison_sm_cropped = comparison_sm(start:ending);

        similarity_framebframe = nan(1, length(comparison_sm_cropped) - d);

        for idx_frame = 1:length(similarity_framebframe)
            v2 = comparison_sm_cropped(idx_frame:idx_frame + d);
            similarity_framebframe(idx_frame) = norm(v1 - v2);
        end

        [closest_similarity, best_frame] = min(similarity_framebframe);
        similarity_value(idx_comparison) = closest_similarity;
        similarity_sort(idx_comparison) = best_frame;
        all_cropped{idx_comparison} = comparison_sm_cropped;
    end

    [sorted_similarity_magnitude, sorted_similarity_idx] = sort(similarity_value);
    best_similarities = sorted_similarity_idx(1:n_selected_comparisons);


    for idx_tops = 1:length(best_similarities)

        comparison_sm_cropped = all_cropped{best_similarities(idx_tops)};

        starting_frame = similarity_sort(best_similarities(idx_tops));
        
        initial_bump(idx_bout, idx_tops) = sum(comparison_sm_cropped(1:starting_frame - 1));

        ending_duration(idx_bout, idx_tops) = length(comparison_sm_cropped(1 : round(y.durations_s(best_similarities(idx_tops)) .* 60))) - length(comparison_sm_cropped(1:starting_frame - 1));

    end

    [init_bump_sorted, init_bump_idx] = sort(initial_bump(idx_bout, :), 'descend');
  
    corr_bout = corrcoef(init_bump_sorted, ending_duration(idx_bout, init_bump_idx) ./60);
    corr_array(idx_bout) = corr_bout(1,2);
end


%%
figure
histogram(corr_array, -1:0.025:1)
apply_generic(gca)

%%
[best_corr, best_corr_idx] = sort(corr_array);

for idx_bout = [best_corr_idx(1:20)]'

    sm = motion_cache(y.fly(idx_bout));
    sm_bout = sm(y.onsets(idx_bout):y.onsets(idx_bout) + round(y.durations_s(idx_bout) .* 60));
    
    offset = length(sm_bout); onset = offset - d;
    v1 = sm_bout(onset:offset);
    similarity_value = nan(1, height(y));
    similarity_sort = nan(1, height(y));

    for idx_comparison = 1:height(y)
        comparison_sm = motion_cache(y.fly(idx_comparison));
        start = y.onsets(idx_comparison);
        ending = y.onsets(idx_comparison) + round(y.durations_s(idx_comparison) .* 60)  + d * percentage_2b_exp;
        comparison_sm_cropped = comparison_sm(start:ending);

        similarity_framebframe = nan(1, length(comparison_sm_cropped) - d);

        for idx_frame = 1:length(similarity_framebframe)
            v2 = comparison_sm_cropped(idx_frame:idx_frame + d);
            similarity_framebframe(idx_frame) = norm(v1 - v2);
        end

        [closest_similarity, best_frame] = min(similarity_framebframe);
        similarity_value(idx_comparison) = closest_similarity;
        similarity_sort(idx_comparison) = best_frame;
        all_cropped{idx_comparison} = comparison_sm_cropped;
    end

    [sorted_similarity_magnitude, sorted_similarity_idx] = sort(similarity_value);
    best_similarities = sorted_similarity_idx(1:n_selected_comparisons);

    fh = figure('color','w','Position',[100, 100, 1600, 800]);
    tiledlayout(2, 3, 'TileSpacing', 'loose', 'Padding', 'loose')
    nexttile(1)
    hold on
    histogram(sorted_similarity_magnitude, 25, 'EdgeColor', 'none' )
    xline(prctile(sorted_similarity_magnitude, 10), 'r--')
    ylabel('Count')
    xlabel('Social Motion')

    apply_generic(gca, 20)
    nexttile(2, [1 2])
    hold on


    for idx_tops = 1:length(best_similarities)

        comparison_sm_cropped = all_cropped{best_similarities(idx_tops)};

        starting_frame = similarity_sort(best_similarities(idx_tops));
        
        initial_bump(idx_bout, idx_tops) = sum(comparison_sm_cropped(1:starting_frame - 1));

        ending_duration(idx_bout, idx_tops) = length(comparison_sm_cropped(1 : round(y.durations_s(best_similarities(idx_tops)) .* 60))) - length(comparison_sm_cropped(1:starting_frame - 1));

    end

    [init_bump_sorted, init_bump_idx] = sort(initial_bump(idx_bout, :), 'descend');

    i = 0;
    
    for idx_tops = init_bump_idx
        
        i = i + 1;

        comparison_sm_cropped = all_cropped{best_similarities(idx_tops)};
        
        starting_frame = similarity_sort(best_similarities(idx_tops));
        
        plot(-starting_frame + 1: - starting_frame + length(comparison_sm_cropped), comparison_sm_cropped, 'Color', col(i, :))
    end

    plot(1:length(v1), v1, 'k--', 'LineWidth', 3)
    xlabel('Time (frames)')
    ylabel('Social Motion')
    apply_generic(gca, 20)

    nexttile
    scatter(init_bump_sorted, y.durations_s(best_similarities(init_bump_idx), :), 40, col, 'filled')
    
    clim([0 200])
    xlabel('Magnitude Initial Social Motion')
    ylabel('Total Freeze Duration (s)')
    yline(y.durations_s(idx_bout), 'k--', 'LineWidth', 2)
    xline(initial_bump(idx_bout, 1), 'k--', 'LineWidth', 2)
    apply_generic(gca,20)
    xlim([-5 205])
    ylim([-0.5 11.5])

    nexttile
    scatter(init_bump_sorted,  ending_duration(idx_bout, init_bump_idx) ./60, 40, col, 'filled')
    
    clim([0 200])
    xlabel('Magnitude Initial Social Motion')
    ylabel('Duration Post-Template (s)')
    yline(ending_duration(idx_bout, 1) ./60, 'k--', 'LineWidth', 2)
    xline(initial_bump(idx_bout, 1) , 'k--', 'LineWidth', 2)
    apply_generic(gca,20)
    ylim([-0.5 11.5])
    xlim([-5 205])

    corr_bout = corrcoef(init_bump_sorted, ending_duration(idx_bout, init_bump_idx) ./60);
    corr_array(idx_bout) = corr_bout(2);
end
