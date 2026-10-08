id_code = 'imm3_mob3_pc4';
paths_out = path_generator('folder', 'experimental_vars/freeze_durations', 'bouts_id', id_code, 'imfirst', false);

col = cmapper();
thresholds = define_thresholds;
thresholds.le_window_fl = [5 40];
thresholds.le_window_sl = [15 50];

bouts = importdata(fullfile(paths_out.dataset, 'bouts.mat'));
bouts = bouts_formatting(bouts, thresholds);

bouts_spontaneous = data_parser_new(bouts, 'period', 'bsl', 'window', 'le', 'type', 'immobility', 'nloom', 10:20);
bouts_le = data_parser_new(bouts, 'period', 'loom', 'window', 'le', 'type', 'immobility', 'nloom', 2:20, 'min_dur', 30, 'exclude_flies', false);
bouts_all = data_parser_new(bouts, 'min_dur', 30, 'exclude_flies', true);

% Process the data: set thresholds for sm, fs and ln. Set minimum duration.
trunc_point = 30;

bouts_proc = data_parser_new(bouts, 'type', 'immobility', 'period', 'loom', 'window', 'le', 'nloom', 2:20, 'min_dur', trunc_point);
threshold = 70;
[bouts_proc, contact_mask] = impose_contact_threshold(bouts_proc, 'threshold', threshold);
bouts_proc.ends = bouts_proc.onsets + bouts_proc.ending_time;
bouts_proc.durations_s = bouts_proc.ending_time ./ 60;
bouts_proc.durations = bouts_proc.ending_time;

sm_freeze_full = extract_sm_from_bouts(bouts_proc, 'type', 'onlyfreeze', 'output_type', 'mat','norm_factor', 10, 'cache', 'motion_cache');
window_dur = 180;
sm_freeze_ili = extract_sm_from_bouts(bouts_proc, 'type', 'onsets', 'output_type', 'mat','norm_factor', 10, 'cache', 'motion_cache', 'window', [0 window_dur]);

bouts_proc.sm = mean(sm_freeze_full, 2, 'omitnan');
bouts_proc.sm_1s = mean(sm_freeze_ili(:,1:60), 2, 'omitnan');
bouts_proc.sm_2s = mean(sm_freeze_ili(:,1:120), 2, 'omitnan');

bouts_proc = bouts_proc(bouts_proc.durations > trunc_point, :);
bouts_proc = bouts_proc(bouts_proc.sm < 1.5,:);
bouts_proc = bouts_proc(bouts_proc.sm_1s < 1.5,:);
type = 'cumulative';
param = {'sm_2s'};
ls = 'fast';
% bouts_proc.ls = 0* ones(height(bouts_proc), 1);
% bouts_proc.sloom = 25 * ones(height(bouts_proc), 1);
fh = fd_distr_withparam_new('bouts', bouts_proc, 'type', type, 'param', param, 'check_quantiles', true,  'export', false, 'paths', paths_out);%;, 'sloom_to_plot', 'fast'); %, 'avg_ss', 'cum_freeze_time', 'avg_fs_1s_norm', 'n_generated_freezes'})

%% save
%exporter(fh, paths_out, sprintf('%s-%s-%s-ls%s.pdf', type, param{1}, period, ls));
