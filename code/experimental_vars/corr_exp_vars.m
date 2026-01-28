paths_out = path_generator('folder', 'descriptive/fd_durs', 'bouts_id', id_code, 'imfirst', false);

id_code = 'imm2_mob2_pc4';
col = cmapper();
thresholds = define_thresholds;
thresholds.le_window_fl = [5 40];
thresholds.le_window_sl = [15 50];

bouts = importdata(fullfile(paths_out.dataset, 'bouts.mat'));
bouts = bouts_formatting(bouts, thresholds);

bouts_le = data_parser_new(bouts, 'period', 'loom', 'window', 'le', 'type', 'immobility', 'nloom', 2:20, 'min_dur', 30);

fh = figure('color','w','Position',[100, 100, 400, 500]);
var1 = 'nloom';
var2 = 'cum_freeze_time';
var1 = 'avg_sm_freeze_norm';

scatter(bouts_le.(var1), bouts_le.(var2))