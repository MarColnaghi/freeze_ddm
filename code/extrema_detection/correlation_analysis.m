%calculate_corrs()


n = 60;
n_selected_comparisons = 120;

% Select Color Map
col = cmapper();
d = 180;
frames_2b_exp = 120;

code4segm = sprintf('d_%d-2bexp_%d', d, frames_2b_exp);

% Specify Paths
run2analyse = 1;
run = sprintf('run%02d', run2analyse);
paths = path_generator('folder', fullfile('/extrema_detection', 'ac_vs_ed', run, code4segm));


fh = figure('color', 'w', 'Position', [100, 100, 600, 300]);
hold on; ax = gca;

for idx_gen_model = {'ac', 'ed'}
    gen_model = idx_gen_model{1};
    load(fullfile(paths.results, sprintf('struct_%s', gen_model)))
    similarities = nan(1, size(s,2));
    correlations = nan(1, size(s,2));

    for idx_freezes = 1:size(s,2)
        temp_s = s(idx_freezes).boutlist;
        similarities(idx_freezes) = temp_s.closest_similarity(n);
        corrmat = corrcoef(temp_s.summed_motion_b4(1:n), temp_s.rt_post_template(1:n));
        correlations(idx_freezes) = corrmat(2);
    end
    if strcmp(gen_model, 'ac')
        histogram(correlations, -1:0.05:1, 'Normalization', 'pdf', 'FaceColor', col.timevarying_sm, 'EdgeColor', 'none')

    elseif strcmp(gen_model, 'ed')
        histogram(correlations, -1:0.05:1, 'Normalization', 'pdf', 'FaceColor', col.extremadetection, 'EdgeColor', 'none');


    end


end

apply_generic(ax)

%% 

for idx_bout = randi(size(s), 1, 5)

    s_temp = s(idx_bout).boutlist;

    fh_ind_bouts = figure('color','w','Position',[100, 100, 1200, 650]);
    tiledlayout(2, 3, 'TileSpacing', 'loose', 'Padding', 'loose')
    nexttile(1)
    hold on
    histogram(s_temp.closest_similarity, 0:0.25:15, 'FaceColor', 'none', 'EdgeColor', 'k', 'DisplayStyle', 'stairs', 'LineWidth', 2, 'HandleVisibility', 'off')
    temp = s_temp.closest_similarity; temp(temp >= s_temp.closest_similarity(n_selected_comparisons)) = nan;
    histogram(temp, 0:0.25:15, 'FaceColor', 'k', 'EdgeColor', 'none', 'DisplayName', sprintf('Top %d Similar', n_selected_comparisons))
    ylabel('Count')
    xlabel('Distance to Template')
    legend('Box', 'off', 'FontSize', 18);
    apply_generic(gca)

    
    nexttile(2, [1 2])
    
    hold on
    fprintf('bout n:%d\n', idx_bout);
    sm = motion_cache(flies(idx_bout));
    freeze_onset = allfr_onsets(idx_bout);
    freeze_duration = allfr_durations(idx_bout);
    freeze_offset = allfr_offsets(idx_bout);

    sm_bout = sm(freeze_onset:freeze_offset);

    template_onset = length(sm_bout) - d + 1;

    v1 = sm_bout(template_onset:end);


    fill([0 d d 0], [-.35 -.35 15.5 15.5], [0.9 0.9 0.9], 'EdgeColor', 'none', 'HandleVisibility', 'off')

    for idx_tops = 1:n_selected_comparisons

        sm = i;
        comparison_sm_cropped ;

        starting_frame = s_temp.best_frame(idx_tops);

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

    exporter(fh_ind_bouts, paths_loop, sprintf('bout_%d.pdf', idx_bout), 'export', export)
end
