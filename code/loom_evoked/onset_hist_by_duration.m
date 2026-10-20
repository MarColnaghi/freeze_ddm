% Histogram of immobility onsets around the loom, one figure per range of freeze
% duration, one tile per loom speed, and one figure with all the ranges side by
% side, a column per range. A real response shows up as a peak that
% sticks out of the flat background of spontaneous onsets, so the figures are
% for judging where the loom-evoked window should sit and how wide it can be
% before it fills with background. The lines are the window under test, set in
% loom_evoked_data. The black outline is the same histogram for the baseline mock
% looms, so what sticks out of it is what the loom added.

[bouts, thresholds, paths] = loom_evoked_data;

% Freeze duration ranges in frames, both ends included. They sort the bouts by
% how long the freeze really was (duration_raw), not by the time it was observed
% for before a contact. Together they cover every bout loom_evoked_bouts keeps,
% which is at least thresholds.min_dur frames
duration_ranges = [0 11; 12 29; 30 59; 60 239; 240 Inf];

export = true;
separate_figures = true;
combined_figure = true;

% The same bouts as the other plots in this folder, contact-censored, for the
% real looms and for the baseline mock looms
bouts_proc = loom_evoked_bouts(bouts, 'min_dur', 3);
bouts_baseline = loom_evoked_bouts(bouts, 'period', 'bsl', 'min_dur', 3);

paths = path_generator('folder', fullfile('loom_evoked','onset_hist', 'CS'));

if separate_figures
    for idx_range = 1:size(duration_ranges, 1)
        le_threshold_plot(bouts_proc, thresholds, duration_ranges(idx_range, :), bouts_baseline, 'durations');

        if export
            exporter(gcf, paths, sprintf('onset_hist_duration_%d_%g.png', duration_ranges(idx_range, 1), duration_ranges(idx_range, 2)))
        end
    end
end

if combined_figure
    % A column per range, each holding the two loom speeds. The legend and the
    % y tick labels are only drawn on the first column, the y axis is shared
    fh = figure('color','w','Position',[100, 100, 1800, 500]);
    tl_all = tiledlayout(fh, 1, size(duration_ranges, 1), 'TileSpacing', 'compact', 'Padding', 'compact');

    for idx_range = 1:size(duration_ranges, 1)
        tl = le_threshold_plot(bouts_proc, thresholds, duration_ranges(idx_range, :), bouts_baseline, 'durations', ...
            'parent', tl_all, 'font_size', 20, 'legend', idx_range == 1, 'yaxis', idx_range == 1 );
        tl.Layout.Tile = idx_range;
    end
    xlabel(tl_all, 'Immobility Onsets (frames)', 'FontSize', 18, 'FontName', 'Helvetica-Light')

    if export
        exporter(fh, paths, 'onset_hist_duration_all.pdf')
    end
end
