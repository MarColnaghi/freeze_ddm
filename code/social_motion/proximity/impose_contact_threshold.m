function bouts_no_contact = impose_contact_threshold(bouts_proc, threshold)

inizio = 1;
fine = 630;

distance_freeze = extract_sm_from_bouts(bouts_proc, 'type', 'onlyfreeze', 'output_type', 'mat', 'cache', 'mindist_cache', 'align', 'onset');
% distance_freeze = extract_sm_from_bouts(bouts_proc, 'type', 'onsets', 'output_type', 'mat', 'cache', 'mindist_cache', 'window', [0 630]);

distance_freeze_cropped = distance_freeze(:, inizio:fine);

contact_mask = any(distance_freeze_cropped <= threshold, 2);
bouts_no_contact = bouts_proc(~contact_mask, :);