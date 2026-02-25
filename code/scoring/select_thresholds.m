
param = {'avg_fs_1s'};
id_code = 'imm2_mob2_pc4';
paths = path_generator('bouts_id', id_code, 'imfirst', false, 'folder', 'descriptive/fd_distr');

thresholds = define_thresholds;
thresholds.le_window_fl = [5 40];
thresholds.le_window_sl = [15 50];

bouts = importdata(fullfile(paths.dataset, 'bouts.mat'));
bouts = bouts_formatting(bouts, thresholds);
bouts_proc = data_parser_new(bouts, 'period', 'loom', 'window', 'le', 'type', 'immobility', 'nloom', 2:20, 'min_dur', 30);
fd_distr_withparam_new('bouts', bouts_proc, 'type', 'ecdf', 'param', param)


thresholds = define_thresholds;
%thresholds.le_window_fl = [10 40];
%thresholds.le_window_sl = [20 50];

bouts = importdata(fullfile(paths.dataset, 'bouts.mat'));
bouts = bouts_formatting(bouts, thresholds);
bouts_proc = data_parser_new(bouts, 'type', 'immobility', 'period', 'loom', 'window', 'le', 'nloom', 2:20, 'min_dur', 30);
fd_distr_withparam_new('bouts', bouts_proc, 'type', 'ecdf', 'param', param)


