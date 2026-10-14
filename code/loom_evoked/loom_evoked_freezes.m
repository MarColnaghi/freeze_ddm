function freezes = loom_evoked_freezes(bouts, varargin)
% The loom-evoked freezes, the bouts the analyses of freeze duration start from.
%
%   freezes = loom_evoked_freezes(bouts)
%   freezes = loom_evoked_freezes(bouts, 'already_frozen', true, 'nloom', 1:20)
%
% They are the bouts of loom_evoked_bouts (immobility of at least min_dur frames
% in the loom period, flies that freeze a lot at baseline excluded,
% contact-censored, seen for at least min_dur frames before any contact) whose
% onset falls in the loom-evoked window of their loom speed, the le column that
% bouts_formatting sets from define_thresholds. durations, durations_s and ends
% hold the time each freeze was observed for, and is_censored says whether it was
% cut short by a contact, by the end of the loom cycle or by the end of the
% recording, which is what ecdf takes together with them. The duration the freeze
% really had is duration_raw.
%
% The window has some background in it. Counting onsets from 12 frames, about one in
% ten would be there without the loom (loom_evoked/split_half_check). Among the
% freezes of at least 30 frames it is about 4%, 2-6% by group, and those are short
% ones: the same selection on the baseline mock looms gives 186 freezes against 4,677
% in the loom period, and only 4-7% of them last more than 3 s.
%
% bouts is the table with bouts_formatting applied, as loom_evoked_data returns it.
%
% Options
%   already_frozen  false  drop the bouts of looms that found the fly already
%                          frozen (frozen_start), as the analyses of freeze
%                          duration have always done. The loom-evoked
%                          diagnostics keep them, so true gives that selection
%   nloom           2:20   loom numbers to keep, the first loom is left out as in
%                          the older analyses
%   min_dur         []     minimum duration in frames, thresholds.min_dur if empty
%   contact_threshold  []  distance to the nearest fly below which a freeze is cut
%                          as a collision, 70 (impose_contact_threshold) if empty

opt = inputParser;
addParameter(opt, 'already_frozen', false, @(x) islogical(x) && isscalar(x));
addParameter(opt, 'nloom', 2:20, @isnumeric);
addParameter(opt, 'min_dur', [], @(x) isempty(x) || (isnumeric(x) && isscalar(x)));
addParameter(opt, 'contact_threshold', [], @(x) isempty(x) || (isnumeric(x) && isscalar(x)));
parse(opt, varargin{:});

freezes = loom_evoked_bouts(bouts, 'min_dur', opt.Results.min_dur, 'contact_threshold', opt.Results.contact_threshold);
freezes = freezes(freezes.le == 1, :);

if ~opt.Results.already_frozen
    freezes = freezes(freezes.frozen_start == 0, :);
end

freezes = freezes(ismember(freezes.nloom, opt.Results.nloom), :);
end
