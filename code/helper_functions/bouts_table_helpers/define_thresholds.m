function thresholds = define_thresholds(varargin)

opt = inputParser;
% Loom-evoked window, onset latency in frames after the loom, for the slow (sl)
% and the fast (fl) loom. These are the original windows, kept on purpose over the
% core band learned from the excess of loom over baseline onsets (the run of
% latencies around the main peak where at least 75% of the onsets are in excess of
% the baseline mock looms), which loom_evoked/split_half_check finds to be [26 43]
% and [14 28], stable to 3 frames across random halves of the flies. At min_dur the
% original windows hold 75% (slow) and 86% (fast) of all evoked onsets, and 0.905
% and 0.885 of the onsets in them are evoked (purity), against 72% and 79% and
% 0.910 and 0.905 for the learned core. The inclusive windows holding 90% of the
% excess are [5 43] and [8 31].
% min_dur is the shortest immobility bout, in frames, that counts as a freeze: the
% shortest one at which the learned core window of both speeds is at least 90% pure.
addParameter(opt, 'le_window', struct('le_window_sl', [26 45], 'le_window_fl', [13 32]));
addParameter(opt, 'min_dur', 12);
parse(opt, varargin{:});

le_window = opt.Results.le_window;

thresholds.pc = 4; thresholds.fill_in_mob = 3; thresholds.fill_in_imm = 3;
thresholds.still_vel = 20; % mm/s, faster frames are moving whatever pixelchange says, as in dataset_resolved

thresholds.le_window_sl = le_window.le_window_sl;
thresholds.le_window_fl = le_window.le_window_fl;
thresholds.min_dur = opt.Results.min_dur;

thresholds.time_window = [120 540];
thresholds.sp_window = [45 75; 100 130; 150 180; 210 240; 300 330; 390 420; -80 -50] - 30;
thresholds.freeze_dur = 1;
thresholds.mob_windows(2,:) = [80 520];