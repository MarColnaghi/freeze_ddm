
type = 'loom';
str = 'immobility';
plot_style = 'bar3';
zoom_flag =  false;
export = false;
bouts_proc = data_parser_new(bouts, 'min_dur', 12);

% Adds censored_contacts / censored_loom / is_censored / ending_time.
% Slow: loads the mindist cache and walks every bout, so keep bouts_proc
% around rather than rebuilding it each time
bouts_proc = impose_contact_threshold(bouts_proc);

hist2d_durs_onsets('bouts', bouts_proc, 'thresholds', thresholds, 'type', type, 'str', str, ...
    'plot_style', plot_style, 'zoom', zoom_flag, 'export', export, 'split_by', 'sloom', ...
    'censoring', 'is_censored')
