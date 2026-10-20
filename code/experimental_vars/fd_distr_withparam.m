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

% Split latency of each loom speed, the dashed line in le_threshold_plot, in frames
% after the loom. late_onset is 1 for the freezes whose onset is after it
speeds = [25 50];
split_latency = [mean(thresholds.le_window_sl) - 1.5, mean(thresholds.le_window_fl) - 2.5];
split_latency_row = reshape(split_latency(1 + (bouts_proc.sloom == 50)), [], 1); % column, indexing a row with a column gives a row
bouts_proc.late_onset = double(bouts_proc.onsets_loomaligned > split_latency_row);

% bouts_proc = bouts_proc(bouts_proc.sm < 1.5,:);
% bouts_proc = bouts_proc(bouts_proc.sm_1s < 1.5,:);
%bouts_proc = bouts_proc(bouts_proc.jump_start == 0, :);
type = 'cumulative';

% Control variables: loom speed, number of moving flies and loom number. Measured
% covariates can be added, 'sm_1s' for the social motion in the first second of the
% freeze and 'fs' for the focal speed before it

params = {'avg_fs_1s'};

fh = fd_distr_withparam_new('bouts', bouts_proc, 'type', type, 'params', params, 'export', false, 'paths', paths_out, 'scale', 'log', 'control_var', true);
fh = fd_distr_withparam_new('bouts', bouts_proc, 'type', type, 'params', params, 'export', false, 'paths', paths_out, 'scale', 'linear', 'control_var', true);

fh = fd_distr_withparam_by_moving_flies('bouts', bouts_proc,'type', type, 'params', params, 'export', false, 'paths', paths_out);

summary = fd_summary_by_control(bouts_proc, params);

%% Early versus late onset within the loom-evoked window
% Does how long the freeze lasts depend on which side of the split latency its
% onset falls, for each loom speed? The split is the dashed line in le_threshold_plot,
% in frames after the loom. Early is onset <= split, late is onset > split. Kaplan-Meier
% curves of the observed duration (contact and the end of the loom cycle as censoring)
% and a Cox model of the hazard of the freeze ending, late against early (HR < 1: late
% onsets last longer), alone and adjusted for loom number and moving flies. The window
% is the same for both groups, but a freeze that starts later has less time left before
% the end of the loom cycle, which the censoring takes care of.

split_summary = table();
for k = 1:numel(speeds)
    T = bouts_proc(bouts_proc.sloom == speeds(k), :);
    late = T.onsets_loomaligned > split_latency(k);
    cens = logical(T.is_censored);

    med = nan(1, 2);
    for g = 0:1
        in_g = late == g;
        [F, x] = ecdf(T.durations_s(in_g), 'Censoring', cens(in_g), 'Function', 'survivor');
        i = find(F <= 0.5, 1);
        if ~isempty(i), med(g + 1) = x(i); end
    end

    [~, ~, ~, st] = coxphfit(double(late), T.durations_s, 'Censoring', cens);
    [~, ~, ~, st_adj] = coxphfit([double(late), T.nloom, T.moving_flies], T.durations_s, 'Censoring', cens);
    row = table(speeds(k), split_latency(k), sum(~late), sum(late), med(1), med(2), exp(st.beta(1)), st.p(1), ...
        exp(st_adj.beta(1)), st_adj.p(1), 'VariableNames', {'sloom', 'split', 'n_early', 'n_late', ...
        'median_early_s', 'median_late_s', 'HR_late', 'p', 'HR_late_adj', 'p_adj'});
    split_summary = [split_summary; row]; %#ok<AGROW>
end
disp(split_summary)

%% save
%exporter(fh, paths_out, sprintf('%s-%s-%s-ls%s.pdf', type, param{1}, period, ls));
