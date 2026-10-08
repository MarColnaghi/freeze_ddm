function bouts_proc = build_freeze_survival(bouts_proc, varargin)
% BUILD_FREEZE_SURVIVAL  Contact-aware survival data and social-motion
% covariates for freeze bouts, in one pass.
%
%   bouts_proc = build_freeze_survival(bouts_proc, 'name', value, ...)
%
% bouts_proc is the table from data_parser_new (selection of which bouts to
% analyse stays in the calling script). This function then
%
%   1. right-censors each bout at its first contact (impose_contact_threshold),
%   2. applies the entry rule: a bout is kept only if it is still frozen and
%      contact-free beyond `landmark` frames, so every kept bout is at risk
%      from the landmark on (left truncation / landmark analysis),
%   3. computes the social-motion covariates over the requested windows, using
%      only frames the bout was actually observed cleanly.
%
% Options
%   contact_threshold  70     distance (px) below which a neighbour counts as
%                             a contact. One value for the whole analysis
%   trunc_point        30     minimum bout length (frames) used upstream in
%                             data_parser_new; the landmark may not be below it
%   landmark           []     frames. Entry time and end of the exposure
%                             window. Default trunc_point
%   windows            {}     N x 2 cell, {name, [first last]}, offsets in
%                             frames from bout onset, inclusive. Negative
%                             offsets are pre-onset. Each gives a column
%                             sm_<name>. Default: 'pre' = [-25 -1] and
%                             'early' = [0 landmark-1]
%   partial            false  false: a window that is not fully observed
%                             (bout or contact ends inside it) gives NaN.
%                             true: average the observed part, which mixes
%                             window lengths across bouts
%   guard              0      frames before the contact (or bout end) also
%                             dropped from windows, to keep the neighbour's
%                             approach out of the covariate
%   sm_max             Inf    drop bouts with sm_<sm_var> above this. This is
%                             selection on the exposure, so report it
%   sm_var             'early'
%   norm_factor        10     divides the raw social motion
%   verbose            true
%
% Columns added
%   t_obs      observed time in frames from onset: the true duration, or the
%              contact frame, or the window edge (631) for a bout that
%              outlasts the loom cycle
%   t_obs_s    t_obs / 60
%   event      true if the freeze ended on its own while observed cleanly,
%              false if right-censored by contact or by the loom cycle
%   sm_<name>  one per window, as described above
% The original columns (durations, ends, ...) are left untouched, unlike the
% older scripts that overwrote them with the censored values. bouts_proc.
% Properties.UserData records the settings, and the contact step is skipped
% when the table already carries the same threshold, so calling this again on
% its own output is cheap.
%
% Why a landmark. A social-motion mean over a window that runs into the
% bout's own risk period conditions on how long the bout lasted: a bout that
% ends at frame 45 gets a 45-frame average and the rest a 60-frame one, and
% the covariate is partly concurrent with the outcome it is meant to predict.
% With the window ending at the landmark, the covariate is fixed before the
% first frame at risk, so it is a clean baseline covariate for Kaplan-Meier
% strata and Cox models, and contact clipping can never touch it.

opt = inputParser;
addParameter(opt, 'contact_threshold', 70);
addParameter(opt, 'trunc_point', 30);
addParameter(opt, 'landmark', []);
addParameter(opt, 'windows', {});
addParameter(opt, 'partial', false);
addParameter(opt, 'guard', 0);
addParameter(opt, 'sm_max', Inf);
addParameter(opt, 'sm_var', 'early');
addParameter(opt, 'norm_factor', 10);
addParameter(opt, 'verbose', true);
parse(opt, varargin{:});

thr         = opt.Results.contact_threshold;
trunc_point = opt.Results.trunc_point;
landmark    = opt.Results.landmark;
if isempty(landmark), landmark = trunc_point; end
if landmark < trunc_point
    error('landmark (%d) is below trunc_point (%d)', landmark, trunc_point);
end
windows = opt.Results.windows;
if isempty(windows)
    windows = {'pre', [-25 -1]; 'early', [0 landmark - 1]};
end
partial = opt.Results.partial;
guard   = opt.Results.guard;

%% Right-censoring at first contact (slow, so skipped when already done)
ud = bouts_proc.Properties.UserData;
done = ismember('ending_time', bouts_proc.Properties.VariableNames) ...
    && isstruct(ud) && isfield(ud, 'contact_threshold') && ud.contact_threshold == thr;
n_in = height(bouts_proc);
if ~done
    bouts_proc = impose_contact_threshold(bouts_proc, 'threshold', thr);
end

bouts_proc.t_obs   = bouts_proc.ending_time;
bouts_proc.t_obs_s = bouts_proc.t_obs ./ 60;
bouts_proc.event   = ~bouts_proc.is_censored;

%% Entry rule: at risk and contact-free beyond the landmark
n_contact_early = sum(bouts_proc.t_obs <= landmark & bouts_proc.censored_contacts);
bouts_proc = bouts_proc(bouts_proc.t_obs > landmark, :);

%% Social motion over each window, restricted to cleanly observed frames
% Frame at offset k after onset is clean if k <= t_obs - 1 - guard: the last
% frame of the bout, or the contact frame itself for a contact-censored bout
last_clean = bouts_proc.t_obs - 1 - guard;
for w = 1:size(windows, 1)
    name = windows{w, 1};
    rng  = windows{w, 2};
    M = extract_sm_from_bouts(bouts_proc, 'type', 'onsets', 'output_type', 'mat', ...
        'window', rng, 'norm_factor', opt.Results.norm_factor, 'cache', 'motion_cache');
    off = rng(1):rng(2);
    clean = off <= last_clean;
    M(~clean) = NaN;
    v = mean(M, 2, 'omitnan');
    if ~partial
        v(~all(clean, 2)) = NaN;
    end
    bouts_proc.(['sm_' name]) = v;
end

%% Trim on the exposure, if asked
n_trim = 0;
if isfinite(opt.Results.sm_max)
    var = ['sm_' opt.Results.sm_var];
    keep = bouts_proc.(var) <= opt.Results.sm_max;   % NaN is dropped too
    n_trim = sum(~keep);
    bouts_proc = bouts_proc(keep, :);
end

bouts_proc.Properties.UserData = struct('contact_threshold', thr, 'trunc_point', trunc_point, ...
    'landmark', landmark, 'windows', {windows}, 'partial', partial, 'guard', guard, ...
    'sm_max', opt.Results.sm_max);

if opt.Results.verbose
    fprintf(['build_freeze_survival: %d bouts in, %d kept (landmark %d f, contact <= %g px).\n' ...
        '  dropped at entry: %d (of which contact within the landmark: %d), by sm_max: %d.\n' ...
        '  kept: %d events, %d contact-censored, %d cycle-censored.\n'], ...
        n_in, height(bouts_proc), landmark, thr, n_in - height(bouts_proc) - n_trim, ...
        n_contact_early, n_trim, sum(bouts_proc.event), sum(bouts_proc.censored_contacts), ...
        sum(bouts_proc.censored_loom));
end
end
