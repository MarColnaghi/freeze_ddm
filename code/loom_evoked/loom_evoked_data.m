function [bouts, thresholds, paths] = loom_evoked_data()
% The bouts table and the loom-evoked window that every script in this folder
% works from, so that changing the window under test is one edit, here.
%
%   [bouts, thresholds, paths] = loom_evoked_data;
%
% bouts has bouts_formatting applied with the window below, so its `le` column
% follows it. Selecting bouts from it is left to each script.
%
% Once a definition is settled, move its windows into define_thresholds and
% leave this on the default.

id_code = 'imm3_mob3_pc4';
paths = path_generator('folder', 'loom_evoked', 'bouts_id', id_code, 'imfirst', false);

% Window under test, onsets in frames after the loom. The default is in
% define_thresholds. Candidates seen across the repo:
%   [5 55] both speeds     [15 50] slow / [5 40] fast     [26 45] slow / [13 32] fast
thresholds = define_thresholds;
% thresholds.le_window_sl = [5 55];
% thresholds.le_window_fl = [5 55];

bouts = importdata(fullfile(paths.dataset, 'bouts.mat'));
bouts = bouts_formatting(bouts, thresholds);
end
