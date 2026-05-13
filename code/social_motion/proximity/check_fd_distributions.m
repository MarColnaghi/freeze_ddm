clearvars

% Load the table first. We will take advantage of an already existing
% dataset.

%% We need to understand what to do with the censoring.

threshold_imm = 2; threshold_mob = 2; threshold_pc = 4; id_code = sprintf('imm%d_mob%d_pc%d', threshold_imm, threshold_mob, threshold_pc);
paths = path_generator('folder', 'social_motion/proximity', 'bouts_id', id_code, 'imfirst', false);
bouts = importdata(fullfile(paths.dataset, 'bouts.mat'));
motion_cache = importdata(fullfile(paths.cache_path, 'motion_cache.mat'));

% Here you can change the window of loom-evoked
thresholds = define_thresholds('le_window', struct('le_window_sl', [5 55], 'le_window_fl', [5 55]));
% thresholds = define_thresholds;

bouts = bouts_formatting(bouts, thresholds);

% Process the data: set thresholds for sm, fs and ln. Set minimum duration.
bouts_proc = data_parser_new(bouts, 'type', 'immobility', 'period', 'loom', 'window', 'le', 'nloom', 2:20, 'min_dur', 30);
threshold = 70;
[bouts_proc, contact_mask] = impose_contact_threshold(bouts_proc, 'threshold', threshold);
bouts_proc.ends = bouts_proc.onsets + bouts_proc.censoring_time;

sm_freeze_full = extract_sm_from_bouts(bouts_proc, 'type', 'onlyfreeze', 'output_type', 'mat','norm_factor', 10, 'cache', 'motion_cache');
bouts_proc.sm = mean(sm_freeze_full, 2, 'omitnan');
bouts_proc = bouts_proc(contact_mask == 0, :);

type = 'cumulative';
param = {'sm'};

fh = fd_distr_withparam_new('bouts', bouts_proc, 'type', type, 'param', param, 'check_quantiles', true,  'export', true, 'paths', paths);
