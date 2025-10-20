%calculate_corrs()


n = 100;

% Select Color Map
col = cmapper(); col.spectral = cbrewer2('Spectral', n);
d = 120;
frames_2b_exp = 60;
extra_frames = d - frames_2b_exp - 1;

code4segm = sprintf('d_%d-2bexp_%d', d, frames_2b_exp);

% Specify Paths
run2analyse = 3;
run = sprintf('run%02d', run2analyse);
paths = path_generator('folder', fullfile('/extrema_detection', 'ac_vs_ed', run));
sim_params = importdata(fullfile(paths.results, 'sim_params.mat'));
bound = sim_params.gt_table.theta_intercept;

fh = figure('color', 'w', 'Position', [100, 100, 600, 300]);
hold on; ax = gca;

for idx_gen_model = {'ac', 'ed'}
    gen_model = idx_gen_model{1};
    load(fullfile(paths.results, code4segm, sprintf('struct_%s', gen_model)))
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

xlabel('Correlation')
apply_generic(ax)
exporter(fh, paths, 'correlations.pdf')

%%
idx_gen_model = {'ac'};
gen_model = idx_gen_model{1};
comparison_sm = importdata(fullfile(paths.results, code4segm, sprintf('comparison_sm_cropped_%s.mat', gen_model)));
load(fullfile(paths.results, code4segm, sprintf('struct_%s', gen_model)))
paths_loop = path_generator('folder', fullfile('/extrema_detection', 'ac_vs_ed', run, code4segm));

for idx_bout = randi(size(s), 1, 10)

    s_temp = s(idx_bout).boutlist;
    s_temp_selected = s_temp(1:n, :);
    s_temp_selected = sortrows(s_temp_selected, 'summed_motion_b4', 'descend');

    fh_ind_bouts = figure('color','w','Position',[100, 100, 1200, 650]);
    tiledlayout(2, 3, 'TileSpacing', 'loose', 'Padding', 'loose')
    nexttile(1)
    hold on
    histogram(s_temp.closest_similarity, 0:0.25:15, 'FaceColor', 'none', 'EdgeColor', 'k', 'DisplayStyle', 'stairs', 'LineWidth', 2, 'HandleVisibility', 'off')
    temp = s_temp.closest_similarity; temp(temp >= s_temp.closest_similarity(n)) = nan;
    histogram(temp, 0:0.25:15, 'FaceColor', 'k', 'EdgeColor', 'none', 'DisplayName', sprintf('Top %d Similar', n))
    ylabel('Count')
    xlabel('Distance to Template')
    legend('Box', 'off', 'FontSize', 18);
    apply_generic(gca)

    nthh = nexttile(2, [1 2]);
    
    hold on
    fprintf('bout n:%d\n', idx_bout);
    sm = s(idx_bout).sm;
    template_onset = length(sm) - d + 1;
    v1 = sm(template_onset:end);
    
    fill([0 d d 0], [-.35 -.35 15.5 15.5], [0.9 0.9 0.9], 'EdgeColor', 'none', 'HandleVisibility', 'off')

    for idx_tops = 1:n
        
        freeze_id = s_temp_selected.idx_freeze(idx_tops);

        comparison_sm_cropped = comparison_sm{freeze_id};
        
        comparison_sm_cropped = comparison_sm_cropped(1:end - extra_frames);

        starting_frame = s_temp_selected.best_frame(idx_tops);

        plot(-starting_frame + 1: - starting_frame + length(comparison_sm_cropped), comparison_sm_cropped, 'Color', [col.spectral(idx_tops, :) 0.8], 'HandleVisibility', 'off')

        scatter(s_temp_selected.rt_post_template(idx_tops) + randn * 0.1, -0.7, 60, [col.spectral(idx_tops, :)], '|', 'LineWidth', 2, 'HandleVisibility', 'off', 'Clipping', 'off', 'MarkerEdgeAlpha', 0.6)
    end

    plot([0 1], [-10 -10], 'Color', col.spectral(1, :), 'LineWidth', 2.5, 'DisplayName', 'Selected Freezes')
    plot(0:length(v1) - 1, v1, 'k-', 'LineWidth', 1, 'DisplayName', 'Template')
    yline(bound, 'k--', 'LineWidth', 1, 'DisplayName', 'Bound')
    xlabel('Template-Aligned Time (frames)')
    ylabel('Social Motion')
    ylim([-.35 5.35]); xlim([-420 420 + d]);
    legend('Box', 'off', 'FontSize', 18);
    apply_generic(gca, 20)
    set(gca ,'Layer', 'Top')
    ax_plots = gcf;
    ax_chil = nthh.Children; 

    nexttile
    hold on
    scatter(s_temp_selected.summed_motion_b4,  s_temp_selected.rt_post_template, 40, col.spectral, 'filled')
    % scatter(initial_bump(idx_bout, 1), ending_duration(idx_bout, 1) ./60, 60, 'k', 'filled')

    clim([0 200])
    xlabel({'Pre-Template', 'Accumulated Soc. Mot.'})
    ylabel('Duration Post-Template (s)')
   % yline(ending_duration(idx_bout, 1) ./60, 'k--', 'LineWidth', 2, 'Label', 'Template', 'FontSize', 18)
   % xline(initial_bump(idx_bout, 1) , 'k--', 'LineWidth', 2)
    apply_generic(gca,20)
    ylim([-0.5 700.5])
    xlim([-50 750])

    nthh = nexttile(5, [1 2]);
    
    hold on
    fprintf('bout n:%d\n', idx_bout);
    sm = s(idx_bout).sm;
    template_onset = length(sm) - d + 1;
    v1 = sm(template_onset:end);
    
    fill([0 d d 0], [-.35 -.35 15.5 15.5], [0.9 0.9 0.9], 'EdgeColor', 'none', 'HandleVisibility', 'off')

    for idx_tops = 1:n
        
        freeze_id = s_temp_selected.idx_freeze(idx_tops);

        comparison_sm_cropped = comparison_sm{freeze_id};
        
        comparison_sm_cropped = comparison_sm_cropped(1:end - extra_frames);

        starting_frame = s_temp_selected.best_frame(idx_tops);

        plot(-starting_frame + 1: - starting_frame + length(comparison_sm_cropped), comparison_sm_cropped, 'Color', [col.spectral(idx_tops, :) 0.8], 'HandleVisibility', 'off' , 'LineWidth', 1.2)

        scatter(s_temp_selected.rt_post_template(idx_tops) + randn * 0.1, -0.7, 60, [col.spectral(idx_tops, :)], '|', 'LineWidth', 2, 'HandleVisibility', 'off', 'Clipping', 'off', 'MarkerEdgeAlpha', 0.6)
    end

    plot([0 1], [-10 -10], 'Color', col.spectral(1, :), 'LineWidth', 2.5, 'DisplayName', 'Selected Freezes')
    plot(0:length(v1) - 1, v1, 'k-', 'LineWidth', 1, 'DisplayName', 'Template')
    yline(bound, 'k--', 'LineWidth', 1, 'DisplayName', 'Bound')
    xlabel('Template-Aligned Time (frames)')
    ylabel('Social Motion')
    ylim([-.35 5.35]); xlim([-420 420 + d]);
    legend('Box', 'off', 'FontSize', 18);
    apply_generic(gca, 20)
    set(gca ,'Layer', 'Top')
    xlim([d - 10, d + 10])
    

    % exporter(fh_ind_bouts, paths, sprintf('bout_%d.pdf', idx_bout))
end
