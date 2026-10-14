function le_threshold_plot(bouts, thresholds, frame_range, baseline, dur_var)
% Histogram of immobility onsets around the loom, one tile per loom speed, with
% the loom-evoked window of `thresholds` drawn on top. frame_range selects the
% freezes by duration: a scalar keeps those of at least that many frames, and
% [min max] keeps the ones between the two, both ends included.
%
% baseline is optional: the same selection for the baseline mock looms, drawn as
% an outline over the bars in the same bins and duration range. Both periods
% have the same flies and 20 looms each, so the counts are directly comparable.
%
% dur_var is the bouts column frame_range applies to, 'durations' by default.
% After impose_contact_threshold that column holds the time a bout was observed
% for, so pass 'duration_raw' to sort bouts by how long the freeze really was

if nargin < 4
    baseline = [];
end

if nargin < 5
    dur_var = 'durations';
end

if isscalar(frame_range)
    frame_range = [frame_range Inf];
end

if isinf(frame_range(2))
    range_label = sprintf('duration >= %d', frame_range(1));
else
    range_label = sprintf('duration %d to %d', frame_range(1), frame_range(2));
end

e.quantiles = 5;
col = cmapper([], 5);
fh = figure('color','w','Position',[100,100, 700, 700]);
tl = tiledlayout(fh, 2, 1, 'TileSpacing', 'tight');

i = 0; values = [];

for idx_slooms = [unique(bouts.sloom)]'
    i = i + 1;
    nexttile
    in_range = bouts.(dur_var) >= frame_range(1) & bouts.(dur_var) <= frame_range(2);
    h_handle = histogram(bouts.onsets_loomwin(in_range & bouts.sloom == idx_slooms), -100.5:1:500.5, 'EdgeColor', 'none', 'FaceColor', col.vars.sloom(3*i, :));
    ax(i) = gca;

    if ~isempty(baseline)
        hold on
        in_range_baseline = baseline.(dur_var) >= frame_range(1) & baseline.(dur_var) <= frame_range(2);
        h_baseline = histogram(baseline.onsets_loomwin(in_range_baseline & baseline.sloom == idx_slooms), -100.5:1:500.5, 'DisplayStyle', 'stairs', 'EdgeColor', 'k', 'LineWidth', 1.2);
        if i == 1
            legend([h_handle h_baseline], {'loom', 'baseline'}, 'Location', 'northeast', 'Box', 'off', 'FontSize', 20, 'AutoUpdate', 'off');
        end
    end

    if idx_slooms == 25
        xline(thresholds.le_window_sl, 'LineWidth', 2)
        % ax(i).XAxis.Visible = 'off';
        xticks(sort([0 thresholds.le_window_sl -50 100]))
        text(-40, 400, range_label, 'HorizontalAlignment', 'left', 'FontSize', 20)

    elseif idx_slooms == 50
        xline(thresholds.le_window_fl, 'LineWidth', 2)
        xticks(sort([0 thresholds.le_window_fl -50 100]))
    end

    apply_generic(ax(i), 'font_size', 24);
    values = [values; h_handle.Values(:)];
end

linkaxes(ax(:))
%ylim([-2 round(max(values), -1) + 50])
ylim([-5 405])
xlim([-50 100])
xlabel('Immobility Onsets (frames)')

