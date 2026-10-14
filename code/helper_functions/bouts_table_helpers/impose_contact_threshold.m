function [bouts_proc, contact_mask, is_below_threshold] = impose_contact_threshold(bouts_proc, varargin)

opt = inputParser;
addParameter(opt, 'threshold', 80);
addParameter(opt, 'type', 'onlyfreeze');
parse(opt, varargin{:});

threshold = opt.Results.threshold;
type = opt.Results.type;

inizio = 1;
fine = 631;

switch type
    case 'onlyfreeze'
        distance_freeze = extract_sm_from_bouts(bouts_proc, 'type', 'onlyfreeze', 'output_type', 'mat', 'cache', 'mindist_cache', 'align', 'onset');
    case 'ili'
        distance_freeze = extract_sm_from_bouts(bouts_proc, 'type', 'onsets', 'output_type', 'mat', 'cache', 'mindist_cache', 'window', [0 fine]);
end

distance_freeze_cropped = distance_freeze(:, inizio:fine);

is_below_threshold = distance_freeze_cropped < threshold;
contact_mask = any(is_below_threshold, 2);
[~, first_idx] = max(is_below_threshold, [], 2);

admin_censor_mask = bouts_proc.durations > fine;
ending_points     = min(bouts_proc.durations, fine);
ending_points(contact_mask) = first_idx(contact_mask);

% A bout that runs into the end of the recording (censored_end, set in
% load_flies_new) has no known end, so it is censored too. It is in neither
% censored_contacts nor censored_loom unless one of those applies as well. A
% bouts.mat built before the column existed has no such bouts
if ismember('censored_end', bouts_proc.Properties.VariableNames)
    rec_censor_mask = logical(bouts_proc.censored_end);
else
    rec_censor_mask = false(height(bouts_proc), 1);
end

% Collision censoring takes priority: if both apply, the collision came first
bouts_proc.is_censored                      = contact_mask | admin_censor_mask | rec_censor_mask;
bouts_proc.censored_contacts                = contact_mask;     
bouts_proc.censored_loom                    = admin_censor_mask & ~contact_mask;
bouts_proc.ending_time                      = ending_points;

bouts_proc.ends = bouts_proc.onsets + bouts_proc.ending_time;
bouts_proc.durations_s = bouts_proc.ending_time ./ 60;
bouts_proc.durations = bouts_proc.ending_time;

window_dur = 630;
sm_freeze_ili = extract_sm_from_bouts(bouts_proc, 'type', 'onsets', 'output_type', 'mat','norm_factor', 10, 'cache', 'motion_cache', 'window', [0 window_dur]);
sm_freeze_full = extract_sm_from_bouts(bouts_proc, 'type', 'onlyfreeze', 'output_type', 'mat','norm_factor', 10, 'cache', 'motion_cache');

bouts_proc.sm = mean(sm_freeze_full, 2, 'omitnan');
bouts_proc.sm_1s = mean(sm_freeze_ili(:,1:60), 2, 'omitnan');
bouts_proc.sm_2s = mean(sm_freeze_ili(:,1:120), 2, 'omitnan');

% bouts_proc = bouts_proc(bouts_proc.sm < 1.5,:);
% bouts_proc = bouts_proc(bouts_proc.sm_1s < 1.5,:);