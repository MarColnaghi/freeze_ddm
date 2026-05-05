clearvars

% Load the table first. We will take advantage of an already existing
% dataset.
threshold_imm = 2; threshold_mob = 2; threshold_pc = 4; id_code = sprintf('imm%d_mob%d_pc%d', threshold_imm, threshold_mob, threshold_pc);
paths = path_generator('folder', 'social_motion/proximity', 'bouts_id', id_code, 'imfirst', false);
bouts = importdata(fullfile(paths.dataset, 'bouts.mat'));
motion_cache = importdata(fullfile(paths.cache_path, 'motion_cache.mat'));

% Here you can change the window of loom-evoked
thresholds = define_thresholds('le_window', struct('le_window_sl', [5 55], 'le_window_fl', [5 55]));
thresholds = define_thresholds;

bouts = bouts_formatting(bouts, thresholds);

% Process the data: set thresholds for sm, fs and ln. Set minimum duration.
bouts_proc = data_parser_new(bouts, 'type', 'immobility', 'period', 'loom', 'window', 'le', 'nloom', 2:20, 'min_dur', 30);
inizio = 1;
fine = 630;

distance_freeze = extract_sm_from_bouts(bouts_proc, 'type', 'onlyfreeze', 'output_type', 'mat', 'cache', 'mindist_cache', 'align', 'onset');
distance_freeze_cropped = distance_freeze(:, inizio:fine);

contact_mask = any(distance_freeze_cropped <= 80, 2);
bouts_no_contact = bouts_proc(~contact_mask, :);

type = 'cumulative';
param = {'sm'};

fh = fd_distr_withparam_new('bouts', bouts_no_contact, 'type', type, 'param', param, 'check_quantiles', true,  'export', false);%;, 'sloom_to_plot', 'fast'); %, 'avg_ss', 'cum_freeze_time', 'avg_fs_1s_norm', 'n_generated_freezes'})
