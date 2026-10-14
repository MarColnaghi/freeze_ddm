function bouts_proc = loom_evoked_bouts(bouts, varargin)
% The bouts every plot in this folder starts from, so that they all show the
% same freezes.
%
%   bouts_proc = loom_evoked_bouts(bouts)                       % real looms
%   bouts_proc = loom_evoked_bouts(bouts, 'period', 'bsl')      % baseline mock looms
%
% Immobility bouts of more than 10 frames, every loom number, flies that
% freeze a lot at baseline excluded. Bouts of flies that were already frozen
% when the loom arrived (frozen_start) are kept, which is why the period is
% picked here and not through data_parser_new, whose 'loom' option drops them.
%
% Contact censoring is applied with impose_contact_threshold (its default
% threshold, 70), and then the length rule once more on the observed time, so only
% bouts seen for more than 10 frames before any contact stay, as in a landmark
% analysis. That adds is_censored, censored_contacts, censored_loom and
% ending_time, and durations / ends now hold the observed time. The duration the
% freeze really had is kept in duration_raw, for plots that sort bouts by it.
%
% bouts is the table as it comes out of loom_evoked_data.
%
% data_parser_new also drops bouts with sm above 1.2 or fs above 2, which is
% part of the shared pipeline and applies here too.

opt = inputParser;
addParameter(opt, 'period', 'loom', @(x) any(strcmp(x, {'loom', 'bsl'})));
parse(opt, varargin{:});

min_dur = 11;      % frames, data_parser_new keeps durations >= min_dur

bouts_proc = data_parser_new(bouts, 'type', 'immobility', 'min_dur', min_dur, 'exclude_flies', true);

switch opt.Results.period
    case 'loom'
        bouts_proc = bouts_proc(bouts_proc.period == 1, :);
    case 'bsl'
        bouts_proc = bouts_proc(bouts_proc.period == 0, :);
end

if isempty(bouts_proc)
    error('loom_evoked_bouts:noBouts', ...
        'No bouts left for period ''%s''. The table passed in has to still contain that period.', opt.Results.period);
end

% Slow, it loads the mindist cache and walks every bout
bouts_proc.duration_raw = bouts_proc.durations;
bouts_proc = impose_contact_threshold(bouts_proc);
bouts_proc = bouts_proc(bouts_proc.durations >= min_dur, :);
end
