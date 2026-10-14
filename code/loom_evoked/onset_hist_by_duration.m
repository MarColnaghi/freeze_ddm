% Histogram of immobility onsets around the loom, one figure per range of freeze
% duration, one tile per loom speed. A real response shows up as a peak that
% sticks out of the flat background of spontaneous onsets, so the figures are
% for judging where the loom-evoked window should sit and how wide it can be
% before it fills with background. The lines are the window under test, set in
% loom_evoked_data. The black outline is the same histogram for the baseline mock
% looms, so what sticks out of it is what the loom added.

[bouts, thresholds, paths] = loom_evoked_data;

% Freeze duration ranges in frames, both ends included. They sort the bouts by
% how long the freeze really was (duration_raw), not by the time it was observed
% for before a contact. Together they cover every bout loom_evoked_bouts keeps,
% which is more than 10 frames
duration_ranges = [11 29; 30 59; 60 119; 120 299; 300 Inf];

export = false;

% The same bouts as the other plots in this folder, contact-censored, for the
% real looms and for the baseline mock looms
bouts_proc = loom_evoked_bouts(bouts);
bouts_baseline = loom_evoked_bouts(bouts, 'period', 'bsl');

for idx_range = 1:size(duration_ranges, 1)
    le_threshold_plot(bouts_proc, thresholds, duration_ranges(idx_range, :), bouts_baseline, 'duration_raw');

    if export
        exporter(gcf, paths, sprintf('onset_hist_duration_%d_%g.png', duration_ranges(idx_range, 1), duration_ranges(idx_range, 2)))
    end
end
