id_code = 'imm3_mob3_pc4';
paths_out = path_generator('folder', 'experimental_vars/freeze_durations', 'bouts_id', id_code, 'imfirst', false);

col = cmapper();
thresholds = define_thresholds;
% thresholds.le_window_fl = [5 40];
% thresholds.le_window_sl = [15 50];

thresholds.le_window_fl = [-10 120];
thresholds.le_window_sl = [-10 120];

bouts = importdata(fullfile(paths_out.dataset, 'bouts.mat'));
bouts = bouts_formatting(bouts, thresholds);

bouts_all = data_parser_new(bouts, 'period', 'loom', 'window', 'all', 'type', 'immobility', 'nloom', 2:20, 'exclude_flies', true, 'min_dur', 30);
ecdf_condt_onset(bouts_all, thresholds)