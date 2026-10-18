% 2D histogram of freeze duration against onset latency around the loom, one
% tile per loom speed. The table and the loom-evoked window come from
% loom_evoked_data, the bouts plotted from loom_evoked_bouts, which already
% applies the contact censoring (slow, it loads the mindist cache), so keep
% bouts_proc around rather than rebuilding it each time.

[bouts, thresholds] = loom_evoked_data;

type = 'loom';
str = 'immobility';
plot_style = 'bar3';
zoom_flag =  false;
export = false;
bouts_proc = loom_evoked_bouts(bouts, 'period', type, 'min_dur', 30);
hist2d_durs_onsets('bouts', bouts_proc, 'thresholds', thresholds, 'type', type, 'str', str, ...
    'plot_style', plot_style, 'zoom', zoom_flag, 'export', true, 'split_by', 'sloom')
