
function correlation_initbump_rt


%% Parse inputs
p = inputParser;
addParameter(p, 'export', false);
addParameter(p, 'type', 'surrogate');
addParameter(p, 'selection', 500);
addParameter(p, 'gen_model', '');  % Plot all trials if true
addParameter(p, 'all', false);  % Plot all trials if true
addParameter(p, 'align_to', 'offset');  % 'offset' (old behavior) or 'onset'


close all

exporting = false;

% Model
model = 'dddm2';
select_run = 'run03';
gen_data = 'fr';

% Specify Paths
paths = path_generator('folder', fullfile('/extrema_detection', model));

% Load the motion ts
sim_params.motion_cache_path = fullfile(paths.dataset, 'motion_cache.mat');
load(sim_params.motion_cache_path)

% Parameters
d = 180; L = d + 1;
percentage_2b_exp = 1 - 1/3; 
n_selected_comparisons = 240;

% Calculations
extra_frames = round(percentage_2b_exp * d) - 1;

% Select Color Map
col = cbrewer2('Spectral', n_selected_comparisons);

fh = figure('color','w','Position',[100, 100, 600, 250]);
hold on

i = 0;

for idx_surrogate_run = 13:16

    i = i + 1;

    clear y
    
    load(sprintf('/Users/marcocolnaghi/PhD/freeze_ddm/model_results/momentary_integration/surrogate/dddm2/run%d/sims_tv/y.mat', idx_surrogate_run))
    paths = path_generator('folder', fullfile('/extrema_detection', model));

    folder_name = fullfile(paths.fig, sprintf('d_%d-2bexp_%.2f-ncomp_%d', d, 1 - percentage_2b_exp, n_selected_comparisons), sprintf('run%d', idx_surrogate_run));
    mkdir(folder_name);
    paths.fig = folder_name;

    % Select only freezes with specific durations
    y = y(y.durations_s > d/60 & y.durations_s <= 10.5, :);

    flies = y.fly;
    allfr_onsets = y.onsets;
    allfr_durations = round(y.durations_s .* 60);
    allfr_offsets = allfr_onsets + allfr_durations;

    % Preallocate once before the outer loop
    initial_bump    = nan(height(y), n_selected_comparisons);
    ending_duration = nan(height(y), n_selected_comparisons);
    ending_duration_sorted = nan(height(y), n_selected_comparisons);
    
    corr_array_surr = nan(height(y), 1);

    for idx_bout = 1:height(y)

        fprintf('running bout #%d out of #%d \n', idx_bout, height(y))

        sm = motion_cache(flies(idx_bout));
        freeze_onset = allfr_onsets(idx_bout);
        freeze_duration = allfr_durations(idx_bout);

        sm_bout = sm(freeze_onset:freeze_onset + freeze_duration);

        template_onset = length(sm_bout) - d + 1;
        
        v1 = sm_bout(template_onset:end);
        
        similarity_value = nan(1, height(y));
        similarity_sort = nan(1, height(y));
        
        all_cropped = cell(1, height(y));

        for idx_comparison = 1:height(y)
            
            comparison_sm = motion_cache(flies(idx_comparison));

            start = allfr_onsets(idx_comparison);
            ending = allfr_offsets(idx_comparison) + extra_frames;
            comparison_sm_cropped = comparison_sm(start:ending);

            similarity_framebframe = nan(1, length(comparison_sm_cropped) - d + 1);

           
            dots = conv(comparison_sm_cropped, flipud(v1), 'valid');
            ssq  = conv(comparison_sm_cropped.^2, ones(d, 1), 'valid');
            v1sq = sum(v1.^2);
           
            similarity_framebframe_vectorized = sqrt(max(v1sq + ssq - 2*dots, 0));

            [closest_similarity, best_frame] = min(similarity_framebframe_vectorized);
            similarity_value(idx_comparison) = closest_similarity;
            similarity_sort(idx_comparison) = best_frame;
            all_cropped{idx_comparison} = comparison_sm_cropped;

        end

        [sorted_similarity_magnitude, best_similarities] = mink(similarity_value, n_selected_comparisons);

        for idx_tops = 1:length(best_similarities)

            comparison_sm_cropped = all_cropped{best_similarities(idx_tops)};

            starting_frame = similarity_sort(best_similarities(idx_tops));

            initial_bump(idx_bout, idx_tops) = sum(comparison_sm_cropped(1:starting_frame - 1));

            ending_duration(idx_bout, idx_tops) = length(comparison_sm_cropped(1 : round(y.durations_s(best_similarities(idx_tops)) .* 60))) - length(comparison_sm_cropped(1:starting_frame - 1));

        end

        [init_bump_sorted, init_bump_idx] = sort(initial_bump(idx_bout, :), 'descend');

        corr_bout = corrcoef(init_bump_sorted, ending_duration(idx_bout, init_bump_idx) ./60);

        ending_duration_sorted(idx_bout, :) =  ending_duration(idx_bout, init_bump_idx);

        corr_array_surr(idx_bout) = corr_bout(1,2);
    end
    
    figure(fh)
    if i == 1
        histogram(corr_array_surr, -1:0.025:1, 'Normalization', 'pdf', 'DisplayStyle', 'stairs', 'DisplayName', 'surrogate')
    else
        histogram(corr_array_surr, -1:0.025:1, 'Normalization', 'pdf',  'DisplayStyle', 'stairs', 'HandleVisibility', 'off')
    end

    imagesc_alltrials(ending_duration_sorted, n_selected_comparison, paths)
end

figure(fh)
apply_generic(gca)

% Here we do it with the empirical

clear y

% Load important files
paths = path_generator('folder', fullfile('/fitting_freezes/le', model, select_run));
load(fullfile(paths.results, 'surrogate.mat'));
y = surrogate;
paths = path_generator('folder', fullfile('/extrema_detection', model));

folder_name = fullfile(paths.fig, sprintf('d_%d-2bexp_%.2f-ncomp_%d', d, 1 - percentage_2b_exp, n_selected_comparisons), 'emp');
mkdir(folder_name);
paths.fig = folder_name;

% Select only freezes with specific durations
y = y(y.durations_s > d/60 & y.durations_s <= 10.5, :);

flies = y.fly;
allfr_onsets = y.onsets;
allfr_durations = round(y.durations_s .* 60);
allfr_offsets = allfr_onsets + allfr_durations;

% Preallocate once before the outer loop
initial_bump    = nan(height(y), n_selected_comparisons);
ending_duration = nan(height(y), n_selected_comparisons);
ending_duration_sorted = nan(height(y), n_selected_comparisons);

corr_array_emp = nan(height(y), 1);

for idx_bout = 1:height(y)

    fprintf('running bout #%d out of #%d \n', idx_bout, height(y))

    sm = motion_cache(flies(idx_bout));
    freeze_onset = allfr_onsets(idx_bout);
    freeze_duration = allfr_durations(idx_bout);

    sm_bout = sm(freeze_onset:freeze_onset + freeze_duration);

    template_onset = length(sm_bout) - d + 1;

    v1 = sm_bout(template_onset:end);

    similarity_value = nan(1, height(y));
    similarity_sort = nan(1, height(y));

    all_cropped = cell(1, height(y));

    for idx_comparison = 1:height(y)

        comparison_sm = motion_cache(flies(idx_comparison));

        start = allfr_onsets(idx_comparison);
        ending = allfr_offsets(idx_comparison) + extra_frames;
        comparison_sm_cropped = comparison_sm(start:ending);

        similarity_framebframe = nan(1, length(comparison_sm_cropped) - d + 1);


        dots = conv(comparison_sm_cropped, flipud(v1), 'valid');
        ssq  = conv(comparison_sm_cropped.^2, ones(d, 1), 'valid');
        v1sq = sum(v1.^2);

        similarity_framebframe_vectorized = sqrt(max(v1sq + ssq - 2*dots, 0));

        [closest_similarity, best_frame] = min(similarity_framebframe_vectorized);
        similarity_value(idx_comparison) = closest_similarity;
        similarity_sort(idx_comparison) = best_frame;
        all_cropped{idx_comparison} = comparison_sm_cropped;

    end

    [sorted_similarity_magnitude, best_similarities] = mink(similarity_value, n_selected_comparisons);

    for idx_tops = 1:length(best_similarities)

        comparison_sm_cropped = all_cropped{best_similarities(idx_tops)};

        starting_frame = similarity_sort(best_similarities(idx_tops));

        initial_bump(idx_bout, idx_tops) = sum(comparison_sm_cropped(1:starting_frame - 1));

        ending_duration(idx_bout, idx_tops) = length(comparison_sm_cropped(1 : round(y.durations_s(best_similarities(idx_tops)) .* 60))) - length(comparison_sm_cropped(1:starting_frame - 1));

    end

    [init_bump_sorted, init_bump_idx] = sort(initial_bump(idx_bout, :), 'descend');

    corr_bout = corrcoef(init_bump_sorted, ending_duration(idx_bout, init_bump_idx) ./60);

    ending_duration_sorted(idx_bout, :) =  ending_duration(idx_bout, init_bump_idx);

    corr_array_emp(idx_bout) = corr_bout(1,2);
end

figure(fh)
histogram(corr_array_emp, -1:0.025:1, 'Normalization', 'pdf', 'EdgeColor', 'none', 'FaceColor', 'k', 'DisplayName', 'empirical')
xlabel('Correlation')
ylabel('density')
lh = legend('Box', 'off', 'FontSize', 20);

paths = path_generator('folder', fullfile('/extrema_detection', model));
exporter(fh, paths, 'correlation_distr.pdf')

imagesc_alltrials(ending_duration_sorted, n_selected_comparison, paths)

folder_name = fullfile(paths.fig, sprintf('d_%d-2bexp_%.2f-ncomp_%d', d, 1 - percentage_2b_exp, n_selected_comparisons), 'emp');
mkdir(folder_name);
paths.fig = folder_name;
[best_corr, best_corr_idx] = sort(corr_array_emp);

for idx_bout = [best_corr_idx(1:5); best_corr_idx(end - 4:end)]'

    fprintf('running bout #%d out of #%d \n', idx_bout, height(y))

    sm = motion_cache(flies(idx_bout));
    freeze_onset = allfr_onsets(idx_bout);
    freeze_duration = allfr_durations(idx_bout);

    sm_bout = sm(freeze_onset:freeze_onset + freeze_duration);

    template_onset = length(sm_bout) - d + 1;

    v1 = sm_bout(template_onset:end);

    similarity_value = nan(1, height(y));
    similarity_sort = nan(1, height(y));

    all_cropped = cell(1, height(y));

    for idx_comparison = 1:height(y)

        comparison_sm = motion_cache(flies(idx_comparison));

        start = allfr_onsets(idx_comparison);
        ending = allfr_offsets(idx_comparison) + extra_frames;
        comparison_sm_cropped = comparison_sm(start:ending);

        similarity_framebframe = nan(1, length(comparison_sm_cropped) - d + 1);


        dots = conv(comparison_sm_cropped, flipud(v1), 'valid');
        ssq  = conv(comparison_sm_cropped.^2, ones(d, 1), 'valid');
        v1sq = sum(v1.^2);

        similarity_framebframe_vectorized = sqrt(max(v1sq + ssq - 2*dots, 0));

        [closest_similarity, best_frame] = min(similarity_framebframe_vectorized);
        similarity_value(idx_comparison) = closest_similarity;
        similarity_sort(idx_comparison) = best_frame;
        all_cropped{idx_comparison} = comparison_sm_cropped;

    end

    [sorted_similarity_magnitude, best_similarities] = mink(similarity_value, n_selected_comparisons);

    fh = figure('color','w','Position',[100, 100, 1200, 650]);
    tiledlayout(2, 3, 'TileSpacing', 'loose', 'Padding', 'loose')
    nexttile(1)
    hold on
    histogram(similarity_value, 50, 'EdgeColor', 'none' )
    xline(sorted_similarity_magnitude(n_selected_comparisons), 'r--')
    ylabel('Count')
    xlabel('Distance to Template')

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
        
        plot(-starting_frame + 1: - starting_frame + length(comparison_sm_cropped), comparison_sm_cropped, 'Color', [col(i, :) 0.8])
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

    exporter(fh, paths, sprintf('bout_%d.pdf', idx_bout))
end

function imagesc_alltrials(ending_duration_sorted, n_selected_comparison, paths)

fh_imgsc = figure('color','w','Position',[100, 100, 300, 700]);
imagesc(ending_duration_sorted, [0 6 * 60])
xlabel({'Comparisons','(Sorted by Inital Bump)'}); 
xticks([1 n_selected_comparisons])
ylabel('Freezes'); yticks([])
cbh = colorbar('Location', 'northoutside', 'LineWidth', 2, 'FontSize', 18);
cbh.Label.String = 'Post-Template Duration(frames)';

apply_generic(gca, 20)
exporter(fh_imgsc, paths, 'imagesc_rt.pdf')
end