function [bouts_proc, contact_mask, is_below_threshold] = impose_contact_threshold(bouts_proc, varargin)

opt = inputParser;
addParameter(opt, 'threshold', 80);
addParameter(opt, 'type', 'onlyfreeze');
parse(opt, varargin{:});

threshold = opt.Results.threshold;
type = opt.Results.type;

inizio = 1;
fine = 630;

switch type
    case 'onlyfreeze'
        distance_freeze = extract_sm_from_bouts(bouts_proc, 'type', 'onlyfreeze', 'output_type', 'mat', 'cache', 'mindist_cache', 'align', 'onset');
    case 'ili'
        distance_freeze = extract_sm_from_bouts(bouts_proc, 'type', 'onsets', 'output_type', 'mat', 'cache', 'mindist_cache', 'window', [0 630]);
end

distance_freeze_cropped = distance_freeze(:, inizio:fine);

is_below_threshold = distance_freeze_cropped < threshold;
[has_contact, first_idx] = max(is_below_threshold, [], 2);

censoring_points = min(bouts_proc.durations, 630);
censoring_points(has_contact) = first_idx(has_contact);

contact_mask = any(is_below_threshold, 2);
bouts_proc.contacts = contact_mask;
bouts_proc.censoring_time = censoring_points;

distance_freeze = extract_sm_from_bouts(bouts_proc, 'type', 'onsets', 'output_type', 'mat', 'cache', 'mindist_cache', 'window', [0 630]);
is_below_threshold = distance_freeze < threshold;
