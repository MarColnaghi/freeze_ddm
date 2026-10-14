% Freeze-exit curves stratified by social motion (Kaplan-Meier, contact right-censored).
%
% Every bout is still frozen and contact-free at `landmark` frames, and the
% social-motion exposure is the mean over [0, landmark), so the covariate is
% fixed before the first frame at risk. See build_freeze_survival and
% fd_km_by_param for the estimand and its caveats.

id_code = 'imm3_mob3_pc4';
paths_out = path_generator('folder', 'experimental_vars/freeze_durations', 'bouts_id', id_code, 'imfirst', false);

thresholds = define_thresholds;
thresholds.le_window_fl = [5 40];
thresholds.le_window_sl = [15 50];

bouts = importdata(fullfile(paths_out.dataset, 'bouts.mat'));
bouts = bouts_formatting(bouts, thresholds);

trunc_point = 30;
bouts_proc = data_parser_new(bouts, 'type', 'immobility', 'period', 'loom', 'window', 'le', ...
    'nloom', 2:20, 'min_dur', trunc_point);

%% Contact censoring, entry rule and social-motion exposure in one call
landmark = 60;      % frames: 30 keeps the most bouts, 60 and 120 give a longer exposure window
bouts_surv = build_freeze_survival(bouts_proc, 'contact_threshold', 70, 'trunc_point', trunc_point, ...
    'landmark', landmark, 'sm_max', 1.5);

%% Figure
[fh, st] = fd_km_by_param('bouts', bouts_surv, 'param', 'sm_early', 'export', false, 'paths', paths_out, ...
    'name', sprintf('km_sm_landmark%d.pdf', landmark));
