% 3 x 4 grid of freeze-duration distributions: social motion (sm_1s) terciles by
% loom speed and focal speed, Kaplan-Meier with contact / loom-cycle censoring.
% sm_1s is the mean over the first second after onset whether or not the freeze
% ends, on purpose (see impose_contact_threshold), so the freezes are those of
% loom_evoked_freezes with no landmark entry rule.

id_code = 'imm3_mob3_pc4';
paths_out = path_generator('folder', 'experimental_vars/freeze_durations', 'bouts_id', id_code, 'imfirst', false);

thresholds = define_thresholds('min_dur', 30);
bouts = importdata(fullfile(paths_out.dataset, 'bouts.mat'));
bouts = bouts_formatting(bouts, thresholds);
freezes = loom_evoked_freezes(bouts, 'min_dur', thresholds.min_dur, 'contact_threshold', 70);
% freezes = freezes(freezes.sm_1s < 1.5, :);   % the km_by_sm trim; off here as in fd_distr_withparam
speeds = [25 50];
split_latency = [mean(thresholds.le_window_sl) - 1.5, mean(thresholds.le_window_fl) - 2.5];
split_latency_row = reshape(split_latency(1 + (freezes.sloom == 50)), [], 1); % column, indexing a row with a column gives a row
freezes.late_onset = double(freezes.onsets_loomaligned > split_latency_row);

% Rows and columns are any variables: names, or {name, n_levels} for quantiles.
% Default: rows = terciles of sm_1s, columns = focal speed (low, high) x loom speed
[fh, st] = fd_km_grid('bouts', freezes, 'rows', {'sm_1s'}, 'cols', {'avg_fs_1s', 'ls'}, ...
    'type', 'hist', 'export', false, 'paths', paths_out, 'name', 'km_grid_hist.pdf');
% Another cut: focal-speed terciles on the rows, loom speed by loom number on the columns
[fh2, ~] = fd_km_grid('bouts', freezes, 'rows', {{'sm', 3}}, 'cols', {'ls', 'avg_fs_1s'}, ...
    'type', 'hist', 'export', false, 'paths', paths_out, 'name', 'km_grid_fs_ls_nloom.pdf');

[fh2, ~] = fd_km_grid('bouts', freezes, 'rows', {{'sm', 3}}, 'cols', {'ls', 'avg_fs_1s'}, ...
    'type', 'cumulative', 'export', false, 'paths', paths_out, 'name', 'km_grid_fs_ls_nloom.pdf');
