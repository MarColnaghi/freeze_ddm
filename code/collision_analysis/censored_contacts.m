
clearvars

% Load the table first. We will take advantage of an already existing
% dataset.
threshold_imm = 2; threshold_mob = 2; threshold_pc = 4; id_code = sprintf('imm%d_mob%d_pc%d', threshold_imm, threshold_mob, threshold_pc);
paths = path_generator('folder', 'social_motion/clustering', 'bouts_id', id_code, 'imfirst', false);
bouts = importdata(fullfile(paths.dataset, 'bouts.mat'));
motion_cache = importdata(fullfile(paths.cache_path, 'motion_cache.mat'));
mindist_cache = importdata(fullfile(paths.cache_path, 'mindist_cache.mat'));

% Here you can change the window of loom-evoked
thresholds = define_thresholds('le_window', struct('le_window_sl', [0 55], 'le_window_fl', [0 55]));
% thresholds = define_thresholds;

bouts = bouts_formatting(bouts, thresholds);
bouts = bouts_proc();

threshold = 90;

% Process the data: set thresholds for sm, fs and ln. Set minimum duration.
bouts_proc = data_parser_new(bouts, 'type', 'immobility', 'period', 'loom', 'window', 'le', 'nloom', 2:20, 'min_dur', 30);
[~, contact_mask] = impose_contact_threshold(bouts_proc, 'threshold', threshold, 'type', 'onlyfreeze');

dist = extract_sm_from_bouts(bouts_proc, 'type', 'onlyfreeze', 'output_type', 'mat', 'cache', 'mindist_cache', 'align', 'offset');
dist_mask = NaN(size(dist)); 
dist_mask(dist > threshold) = 1;

sm_onlyfreeze = extract_sm_from_bouts(bouts_proc, 'type', 'onlyfreeze', 'output_type', 'mat', 'norm_factor', 10, 'align', 'offset');
sm_only_freeze_withdistancemask = sm_onlyfreeze .* dist_mask;

figure
hold on

contact = median(sm_only_freeze_withdistancemask(contact_mask, :), 'omitnan');
no_contact = median(sm_only_freeze_withdistancemask(~contact_mask, :), 'omitnan');

plot(-numel(contact):1:-1, contact, 'k-')
plot(-numel(no_contact):1:-1, no_contact, 'r-')