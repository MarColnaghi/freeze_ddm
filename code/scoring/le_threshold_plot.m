function le_threshold_plot(bouts, thresholds, frame_threshold, paths, export)

e.quantiles = 5;
col = cmapper([], e);
fh = figure('color','w','Position',[100,100, 700, 700]);
tl = tiledlayout(fh, 2, 1, 'TileSpacing', 'tight');

i = 0; values = [];

for idx_slooms = [unique(bouts.sloom)]'
    i = i + 1;
    nexttile
    h_handle = histogram(bouts.onsets_loomwin(bouts.durations >= frame_threshold  & bouts.sloom == idx_slooms), -100.5:1:500.5, 'EdgeColor', 'none', 'FaceColor', col.vars.sloom(3*i, :));
    ax(i) = gca;

    if idx_slooms == 25
        xline(thresholds.le_window_sl, 'LineWidth', 2)
        % ax(i).XAxis.Visible = 'off';
        xticks(sort([0 thresholds.le_window_sl -50 100]))
        text(-40, 400, sprintf('threshold: duration > %d',frame_threshold), 'HorizontalAlignment', 'left', 'FontSize', 20)

    elseif idx_slooms == 50
        xline(thresholds.le_window_fl, 'LineWidth', 2)
        xticks(sort([0 thresholds.le_window_fl -50 100]))
    end

    apply_generic(ax(i), 24);
    values = [values; h_handle.Values(:)];
end

linkaxes(ax(:))
%ylim([-2 round(max(values), -1) + 50])
ylim([-5 405])
xlim([-50 100])
xlabel('Immobility Onsets (frames)')

if export
    figure_title = sprintf('le_thresholds');

    exporter(fh, paths, [figure_title, '.png'])
    exporter(fh, paths, [figure_title, '.pdf'])
end