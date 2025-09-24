
function correlation_initbump_rt_ed_ac(varargin)


%% Parse inputs
p = inputParser;
addParameter(p, 'export', false);
addParameter(p, 'template_length', '');
addParameter(p, 'frames_2b_exp', '');
addParameter(p, 'n_selected_comparisons', '');  % Plot all trials if true
addParameter(p, 'surrogate_runs', 11:16);  % 'offset' (old behavior) or 'onset'

parse(p, varargin{:});

export = p.Results.export;
d = p.Results.template_length;
frames_2b_exp = p.Results.frames_2b_exp;
n_selected_comparisons = p.Results.n_selected_comparisons;
surrogate_runs = p.Results.surrogate_runs;

% Model: This allows to get the latest run of the fitting procedure.
model = 'dddm2';
select_run = 'run03';
load(fullfile('/Users/marcocolnaghi/PhD/freeze_ddm/model_results/fitting_freezes/le', model, select_run, sprintf('fit_results_%s.mat', model)));

% Specify Paths
paths = path_generator('folder', fullfile('/extrema_detection', model));

% Load the Motion Cache
sim_params.motion_cache_path = fullfile(paths.dataset, 'motion_cache.mat');
load(sim_params.motion_cache_path)

% Calculations
extra_frames = d - frames_2b_exp - 1;

% Select Color Map
col = cbrewer2('Spectral', n_selected_comparisons);

% Already Creates the figure for the histogram
fh = figure('color', 'w', 'Position', [100, 100, 600, 250]);
hold on

i = 0;
surrogate_runs = num2cell(surrogate_runs, 1);
surrogate_runs{end + 1} = 'emp';

for idx_surrogate_run = surrogate_runs

    %fprintf('running surrogate run #%d out of #%d \n', idx_surrogate_run, surrogate_runs(end))
    idx_surrogate_loop = idx_surrogate_run{1};

    i = i + 1;

    if ~ischar(idx_surrogate_loop)

        y_surrogate = importdata(sprintf('/Users/marcocolnaghi/PhD/freeze_ddm/model_results/momentary_integration/surrogate/%s/run%d/sims_tv/y.mat', model, idx_surrogate_loop));
        folder_name = fullfile(sprintf('d_%d-2bexp_%d-ncomp_%d', d, frames_2b_exp, n_selected_comparisons), sprintf('run%d', idx_surrogate_loop));

    else

        paths = path_generator('folder', fullfile('/fitting_freezes/le', model, select_run));
        y_surrogate = importdata(fullfile(paths.results, 'surrogate.mat'));
        folder_name = fullfile(sprintf('d_%d-2bexp_%d-ncomp_%d', d, frames_2b_exp, n_selected_comparisons), 'emp');

    end

    paths_loop = path_generator('folder', fullfile('/extrema_detection', model, folder_name));
    mkdir(paths_loop.fig)

    % Select only freezes with specific durations
    y_surrogate = y_surrogate(y_surrogate.durations_s > d/60 & y_surrogate.durations_s < points.censoring, :);

    flies = y_surrogate.fly;
    allfr_onsets = y_surrogate.onsets;
    allfr_durations = round(y_surrogate.durations_s .* 60);
    allfr_offsets = allfr_onsets + allfr_durations;

    % Preallocate once before the outer loop
    initial_bump    = nan(height(y_surrogate), n_selected_comparisons);
    ending_duration = nan(height(y_surrogate), n_selected_comparisons);
    ending_duration_sorted = nan(height(y_surrogate), n_selected_comparisons);
    
    corr_array_surr = nan(height(y_surrogate), 1);

    for idx_bout = 1:height(y_surrogate)

        sm = motion_cache(flies(idx_bout));
        freeze_onset = allfr_onsets(idx_bout);
        freeze_duration = allfr_durations(idx_bout);

        sm_bout = sm(freeze_onset:freeze_onset + freeze_duration);

        template_onset = length(sm_bout) - d + 1;
        
        v1 = sm_bout(template_onset:end);
        
        similarity_value = nan(1, height(y_surrogate));
        similarity_sort = nan(1, height(y_surrogate));
        
        all_cropped = cell(1, height(y_surrogate));

        for idx_comparison = 1:height(y_surrogate)
            
            comparison_sm = motion_cache(flies(idx_comparison));

            start = allfr_onsets(idx_comparison);
            ending = allfr_offsets(idx_comparison) + extra_frames;
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

            ending_duration(idx_bout, idx_tops) = length(comparison_sm_cropped(1 : round(y_surrogate.durations_s(best_similarities(idx_tops)) .* 60))) - length(comparison_sm_cropped(1:starting_frame - 1));

        end

        [init_bump_sorted, init_bump_idx] = sort(initial_bump(idx_bout, :), 'descend');

        corr_bout = corrcoef(init_bump_sorted, ending_duration(idx_bout, init_bump_idx) ./60);

        ending_duration_sorted(idx_bout, :) =  ending_duration(idx_bout, init_bump_idx);

        corr_array_surr(idx_bout) = corr_bout(1,2);
    end


    figure(fh)

    if i == 1
        histogram(corr_array_surr, -1:0.025:1, 'Normalization', 'pdf', 'DisplayStyle', 'stairs', 'DisplayName', 'surrogate', 'LineWidth', 2)
    else
        if ~ischar(idx_surrogate_loop)
            histogram(corr_array_surr, -1:0.025:1, 'Normalization', 'pdf',  'DisplayStyle', 'stairs', 'HandleVisibility', 'off', 'LineWidth', 2)
        else
            histogram(corr_array_surr, -1:0.025:1, 'Normalization', 'pdf',  'DisplayStyle', 'bar', 'DisplayName', 'emprical', 'FaceColor', 'k', 'EdgeColor', 'none')

        end
    end

    imagesc_alltrials(ending_duration_sorted, n_selected_comparisons, paths_loop)

    [~, best_corr_idx] = sort(corr_array_surr);

    for idx_bout = [best_corr_idx(1:5); best_corr_idx(end - 4:end)]'

    sm = motion_cache(flies(idx_bout));
    freeze_onset = allfr_onsets(idx_bout);
    freeze_duration = allfr_durations(idx_bout);

    sm_bout = sm(freeze_onset:freeze_onset + freeze_duration);

    template_onset = length(sm_bout) - d + 1;

    v1 = sm_bout(template_onset:end);

    similarity_value = nan(1, height(y_surrogate));
    similarity_sort = nan(1, height(y_surrogate));

    all_cropped = cell(1, height(y_surrogate));

    for idx_comparison = 1:height(y_surrogate)

        comparison_sm = motion_cache(flies(idx_comparison));

        start = allfr_onsets(idx_comparison);
        ending = allfr_offsets(idx_comparison) + extra_frames;
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

    [sorted_similarity_magnitude, best_similarities] = mink(similarity_value, n_selected_comparisons);

    fh_ind_bouts = figure('color','w','Position',[100, 100, 1200, 650]);
    tiledlayout(2, 3, 'TileSpacing', 'loose', 'Padding', 'loose')
    nexttile(1)
    hold on
    histogram(similarity_value, 0:0.25:15, 'FaceColor', 'none', 'EdgeColor', 'k', 'DisplayStyle', 'stairs', 'LineWidth', 2, 'HandleVisibility', 'off')
    temp = similarity_value; temp(temp >= sorted_similarity_magnitude(n_selected_comparisons)) = nan;
    histogram(temp, 0:0.25:15, 'FaceColor', 'k', 'EdgeColor', 'none', 'DisplayName', sprintf('Top %d Similar', n_selected_comparisons))

    ylabel('Count')
    xlabel('Distance to Template')
    legend('Box', 'off', 'FontSize', 18);

    apply_generic(gca, 20)
    nexttile(2, [1 2])
    hold on
    fill([0 d d 0], [-.35 -.35 15.5 15.5], [0.9 0.9 0.9], 'EdgeColor', 'none', 'HandleVisibility', 'off')

    for idx_tops = 1:length(best_similarities)

        comparison_sm_cropped = all_cropped{best_similarities(idx_tops)};

        starting_frame = similarity_sort(best_similarities(idx_tops));

        initial_bump(idx_bout, idx_tops) = sum(comparison_sm_cropped(1:starting_frame - 1));

        ending_duration(idx_bout, idx_tops) = length(comparison_sm_cropped(1 : round(y_surrogate.durations_s(best_similarities(idx_tops)) .* 60))) - length(comparison_sm_cropped(1:starting_frame - 1));

    end

    [init_bump_sorted, init_bump_idx] = sort(initial_bump(idx_bout, :), 'descend');

    i = 0;

    for idx_tops = init_bump_idx

        i = i + 1;

        comparison_sm_cropped = all_cropped{best_similarities(idx_tops)};

        starting_frame = similarity_sort(best_similarities(idx_tops));

        plot(-starting_frame + 1: - starting_frame + length(comparison_sm_cropped), comparison_sm_cropped, 'Color', [col(i, :) 0.8], 'HandleVisibility', 'off')
        
        scatter(ending_duration(idx_bout, idx_tops), -0.7, 60, [col(i, :)], '|', 'LineWidth', 2, 'HandleVisibility', 'off', 'Clipping', 'off', 'MarkerEdgeAlpha', 0.6)
    end
    
    plot([0 1], [-10 -10], 'Color', col(1, :), 'LineWidth', 2.5, 'DisplayName', 'Selected Freezes')
    plot(1:length(v1), v1, 'k-', 'LineWidth', 2.5, 'DisplayName', 'Template')
    xlabel('Template-Aligned Time (frames)')
    ylabel('Social Motion')
    ylim([-.35 15.5]); xlim([-420 420 + d]);
    legend('Box', 'off', 'FontSize', 18);
    set(gca ,'Layer', 'Top')

    apply_generic(gca, 20)

    nexttile
    hold on
    scatter(init_bump_sorted, y_surrogate.durations_s(best_similarities(init_bump_idx), :), 40, col, 'filled')
    scatter(initial_bump(idx_bout, 1), y_surrogate.durations_s(best_similarities(1), :), 60, 'k', 'filled')

    clim([0 200])
    xlabel({'Pre-Template', 'Accumulated Soc. Mot.'})
    ylabel('Total Freeze Duration (s)')
    yline(y_surrogate.durations_s(idx_bout), 'k--', 'LineWidth', 2, 'Label', 'Template', 'FontSize', 18)
    xline(initial_bump(idx_bout, 1), 'k--', 'LineWidth', 2)
    apply_generic(gca,20)
    xlim([-5 305])
    ylim([-0.5 11.5])

    nexttile
    hold on
    scatter(init_bump_sorted,  ending_duration(idx_bout, init_bump_idx) ./60, 40, col, 'filled')
    scatter(initial_bump(idx_bout, 1), ending_duration(idx_bout, 1) ./60, 60, 'k', 'filled')

    clim([0 200])
    xlabel({'Pre-Template', 'Accumulated Soc. Mot.'})
    ylabel('Duration Post-Template (s)')
    yline(ending_duration(idx_bout, 1) ./60, 'k--', 'LineWidth', 2, 'Label', 'Template', 'FontSize', 18)
    xline(initial_bump(idx_bout, 1) , 'k--', 'LineWidth', 2)
    apply_generic(gca,20)
    ylim([-0.5 11.5])
    xlim([-5 305])

    exporter(fh_ind_bouts, paths_loop, sprintf('bout_%d.pdf', idx_bout))
   
    end
end

figure(fh)
xlabel('Correlation')
ylabel('density')
legend('Box', 'off', 'FontSize', 18);
xlim([-1 1])
apply_generic(gca, 20)

paths = path_generator('folder', fullfile('/extrema_detection', model));
exporter(fh, paths, 'correlation.pdf')

end



% 
% if idx_bout == 2
%     fcomp = figure('color', 'w', 'Position', [100, 100, 1300, 400]);
%     hold on
%     plot(1:length(v1), v1, 'k', 'LineWidth', 2)
%     plot(length(v1) - freeze_duration:length(v1), sm_bout, 'Color', 'red', 'DisplayName', 'Reference Bout')
%     ax = gca;
%     apply_generic(gca, 20)
%     xticks([length(v1) - freeze_duration; 0; d]);
%     lbl = {sprintf('Freeze Onset'), sprintf('Template Onset'), sprintf('Template End')};
%     ylim([-0.3 1.3]); ax.YAxis.Visible ='off';
%     xticklabels(lbl)
%     ax.XAxis.TickLabelInterpreter = 'tex';
% 
%     plot(comparison_sm_cropped + .5 , 'Color', 'blue')
% end

function imagesc_alltrials(ending_duration_sorted, n_selected_comparisons, paths)

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