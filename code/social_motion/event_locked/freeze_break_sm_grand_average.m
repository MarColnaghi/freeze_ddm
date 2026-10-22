% Social motion aligned to the moment the focal fly breaks its freeze, overall and
% split by the control variables (loom speed, number of moving flies, ...), one
% curve per level with a SEM band over flies. The work is in
% freeze_break_sm_by_param.
%
% The freezes are those of loom_evoked_freezes, as in
% experimental_vars/fd_distr_withparam: loom-evoked, at least thresholds.min_dur
% frames, contact-censored. Only the ones that ended by themselves have a break
% to align to, so the figure is conditional on that (no contact, shorter than the
% observation window): the contact censoring is driven by social motion.

id_code = 'imm3_mob3_pc4';
paths_out = path_generator('folder', 'social_motion/event_locked', 'bouts_id', id_code, 'imfirst', false);

thresholds = define_thresholds('min_dur', 30);

bouts = importdata(fullfile(paths_out.dataset, 'bouts.mat'));
bouts = bouts_formatting(bouts, thresholds);

contact_threshold = 70;
bouts_proc = loom_evoked_freezes(bouts, 'already_frozen', false, 'min_dur', thresholds.min_dur, 'contact_threshold', contact_threshold);

% 'all' is the single grand average. Measured covariates can be added, 'sm_1s'
% for the social motion in the first second of the freeze and 'fs' for the focal speed before it
params = {'moving_flies'};

fh = freeze_break_sm_by_param('bouts', bouts_proc, 'params', params, 'export', true, 'paths', paths_out, 'split_sloom', false);
