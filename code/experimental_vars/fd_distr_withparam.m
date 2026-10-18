% Loom-evoked freeze duration against the control variables: one curve of how long
% the freezes last for each level of each variable, with contact and the end of the
% loom cycle as censoring, and the same in numbers (fd_summary_by_control).
%
% The freezes are those of loom_evoked_freezes: immobility bouts of at least
% thresholds.min_dur frames whose onset falls in the loom-evoked window of
% define_thresholds, contact-censored. The curves are conditional on the freeze
% having lasted the minimum duration.

id_code = 'imm3_mob3_pc4';
paths_out = path_generator('folder', 'experimental_vars/freeze_durations', 'bouts_id', id_code, 'imfirst', false);

col = cmapper();
thresholds = define_thresholds('min_dur', 30);

bouts = importdata(fullfile(paths_out.dataset, 'bouts.mat'));
bouts = bouts_formatting(bouts, thresholds);

% By default the first loom and the bouts of looms that found the fly already
% frozen are left out, as in the analyses before the new definition. 'nloom', 1:20
% and 'already_frozen', true keep them, as the loom-evoked diagnostics do.
%
% The collisions are applied in this call (loom_evoked_freezes > loom_evoked_bouts >
% impose_contact_threshold): a freeze is cut at the first frame the nearest fly is
% closer than contact_threshold, and is_censored says which freezes were cut, by a
% contact or by the end of the loom cycle. The length rule is then applied to the
% time observed before the contact, so a freeze that is touched within min_dur
% frames is not in the analysis. fd_distr_withparam_new and fd_summary_by_control
% take is_censored as the censoring of the curves

contact_threshold = 70;
bouts_proc = loom_evoked_freezes(bouts, 'already_frozen', false, 'min_dur', thresholds.min_dur, 'contact_threshold', contact_threshold);

n_freezes = height(bouts_proc);
n_end_of_record = sum(bouts_proc.is_censored & ~bouts_proc.censored_contacts & ~bouts_proc.censored_loom);
fprintf(['%d loom-evoked freezes, at least %d frames, contact threshold %g: %.1f%% cut by a contact, ' ...
    '%.1f%% by the end of the loom cycle, %.1f%% by the end of the recording\n'], n_freezes, thresholds.min_dur, ...
    contact_threshold, 100 * mean(bouts_proc.censored_contacts), 100 * mean(bouts_proc.censored_loom), 100 * n_end_of_record / n_freezes);

% bouts_proc = bouts_proc(bouts_proc.sm < 1.5,:);
% bouts_proc = bouts_proc(bouts_proc.sm_1s < 1.5,:);
%bouts_proc = bouts_proc(bouts_proc.jump_start == 0, :);
type = 'kde';

% Control variables: loom speed, number of moving flies and loom number. Measured
% covariates can be added, 'sm_1s' for the social motion in the first second of the
% freeze and 'fs' for the focal speed before it

params = {'ls'};

fh = fd_distr_withparam_new('bouts', bouts_proc, 'type', type, 'params', params, 'export', true, 'paths', paths_out, 'scale', 'log', 'control_var', false);
fh = fd_distr_withparam_new('bouts', bouts_proc, 'type', type, 'params', params, 'export', true, 'paths', paths_out, 'scale', 'linear', 'control_var', false);

fh = fd_distr_withparam_by_moving_flies('bouts', bouts_proc,'type', type, 'params', params, 'export', true, 'paths', paths_out);

summary = fd_summary_by_control(bouts_proc, params);

%% save
%exporter(fh, paths_out, sprintf('%s-%s-%s-ls%s.pdf', type, param{1}, period, ls));
