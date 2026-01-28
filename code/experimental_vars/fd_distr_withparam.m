id_code = 'imm3_mob3_pc4';
paths_out = path_generator('folder', 'descriptive/fd_durs', 'bouts_id', id_code, 'imfirst', false);

col = cmapper();
thresholds = define_thresholds;
thresholds.le_window_fl = [5 40];
thresholds.le_window_sl = [15 50];

bouts = importdata(fullfile(paths_out.dataset, 'bouts.mat'));
bouts = bouts_formatting(bouts, thresholds);

paths_out = path_generator('folder', '/spontaneous_process', 'bouts_id', id_code, 'imfirst', false);
bouts_spontaneous = data_parser_new(bouts, 'period', 'bsl', 'window', 'le', 'type', 'immobility', 'nloom', 10:20);
bouts_le = data_parser_new(bouts, 'period', 'loom', 'window', 'le', 'type', 'immobility', 'nloom', 2:20, 'min_dur', 30);
bouts_all = data_parser_new(bouts, 'min_dur', 30);

fd_distr_withparam_new('bouts', bouts_le, 'type', 'ecdf', 'param', {'avg_history_dur'}); %, 'avg_ss', 'cum_freeze_time', 'avg_fs_1s_norm', 'n_generated_freezes'})