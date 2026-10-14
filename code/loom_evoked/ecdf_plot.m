% ECDF of freeze durations, one curve per onset latency, one column of tiles per
% loom speed: the real looms on top and the baseline mock looms below, on the
% same axes and colours. The colour bar greys out the latencies outside the
% loom-evoked window under test (set in loom_evoked_data). The bouts plotted are
% the ones picked by loom_evoked_bouts.

[bouts, thresholds] = loom_evoked_data;

% Frames of onset latency pooled into one curve, 1 is a curve per latency
bin_size = 3;

ecdf_condt_onset(bouts, thresholds, 'bin_size', bin_size)
