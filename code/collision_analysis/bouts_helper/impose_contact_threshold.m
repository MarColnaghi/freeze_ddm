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
        distance_freeze = extract_sm_from_bouts(bouts_proc, 'type', 'onsets', 'output_type', 'mat', 'cache', 'mindist_cache', 'window', [0 fine]);
end

distance_freeze_cropped = distance_freeze(:, inizio:fine);

is_below_threshold = distance_freeze_cropped < threshold;
contact_mask = any(is_below_threshold, 2);
[~, first_idx] = max(is_below_threshold, [], 2);

admin_censor_mask = bouts_proc.durations > fine;   % bout exceeded observation window
ending_points     = min(bouts_proc.durations, fine);
ending_points(contact_mask) = first_idx(contact_mask);

% Collision censoring takes priority: if both apply, the collision came first
bouts_proc.is_censored        = contact_mask | admin_censor_mask;
bouts_proc.contacts           = contact_mask;        % alias kept for compatibility
bouts_proc.censored_collision = contact_mask;
bouts_proc.censored_admin     = admin_censor_mask & ~contact_mask;
bouts_proc.ending_time        = ending_points;
bouts_proc.censoring_time     = ending_points;       % alias kept for compatibility
