function [bouts, soc_mot] = bouts_formatting(bouts, thresholds, soc_mot)

bouts.id = [1:height(bouts)]'; 
soc_mot.id = [1:height(soc_mot)]';

bouts.onsets_loomwin(bouts.onsets_loomaligned <= thresholds.time_window(2)) = bouts.onsets_loomaligned(bouts.onsets_loomaligned <= thresholds.time_window(2));
bouts.onsets_loomwin(bouts.onsets_loomwin <= -thresholds.time_window(1)) = nan;
bouts.nloom_loomwin(bouts.onsets_loomwin < 0) = bouts.nloom(bouts.onsets_loomwin < 0) + 1;
bouts.le = logical(bouts.sloom == 25 & bouts.onsets_loomaligned >= thresholds.le_window_sl(1) & bouts.onsets_loomaligned <= thresholds.le_window_sl(2) | ...
    bouts.sloom == 50 & bouts.onsets_loomaligned >= thresholds.le_window_fl(1) & bouts.onsets_loomaligned <= thresholds.le_window_fl(2));
