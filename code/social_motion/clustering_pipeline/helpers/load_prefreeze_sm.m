function out = load_prefreeze_sm(varargin)
% LOAD_PREFREEZE_SM  Load loom-evoked freeze bouts and the pre-freeze social
% motion matrix, with collision-contaminated frames masked out. Shared setup
% for the clustering_pipeline scripts.
%
%   out = load_prefreeze_sm()                uses the defaults below
%   out = load_prefreeze_sm('name', value)   overrides individual options
%
% Options (defaults match clustering_sm_timeseries.m):
%   le_window         [0 55]   loom-evoked window (sl & fl)
%   nloom             2:20
%   min_dur           30       minimum freeze duration (frames)
%   contact_threshold 0        ili contact-distance threshold (px)
%   window            [0 630]  pre-freeze SM window (frames, from onset)
%   norm_factor       10
%
% Fields of OUT:
%   bouts_proc          processed bout table (contact-free loom-evoked freezes)
%   sm_freeze_full      trials x time pre-freeze social motion (contact-masked)
%   is_below_threshold  trials x time logical, frames in contact
%   inizio, fine        window bounds (frames) used for sm_freeze_full

opt = inputParser;
addParameter(opt, 'le_window', [0 55]);
addParameter(opt, 'nloom', 2:20);
addParameter(opt, 'min_dur', 30);
addParameter(opt, 'contact_threshold', 0);
addParameter(opt, 'window', [0 630]);
addParameter(opt, 'norm_factor', 10);
parse(opt, varargin{:});

% Load the table
threshold_imm = 2; threshold_mob = 2; threshold_pc = 4;
id_code = sprintf('imm%d_mob%d_pc%d', threshold_imm, threshold_mob, threshold_pc);
paths = path_generator('folder', 'social_motion/clustering', 'bouts_id', id_code, 'imfirst', false);
bouts = importdata(fullfile(paths.dataset, 'bouts.mat'));

w = opt.Results.le_window;
thresholds = define_thresholds('le_window', struct('le_window_sl', w, 'le_window_fl', w));
bouts = bouts_formatting(bouts, thresholds);

% Loom-evoked freezes only, drop collision-contaminated bouts
bouts_proc = data_parser_new(bouts, 'type', 'immobility', 'period', 'loom', ...
    'window', 'le', 'nloom', opt.Results.nloom, 'min_dur', opt.Results.min_dur);
[bouts_proc, contact_mask, is_below_threshold] = impose_contact_threshold(bouts_proc, ...
    'threshold', opt.Results.contact_threshold, 'type', 'ili');
bouts_proc         = bouts_proc(contact_mask == 0, :);
is_below_threshold = is_below_threshold(contact_mask == 0, :);

% Pre-freeze social motion (contact frames zeroed out)
inizio = opt.Results.window(1);
fine   = opt.Results.window(2);
sm_freeze_full = extract_sm_from_bouts(bouts_proc, 'type', 'onsets', 'output_type', 'mat', ...
    'window', [inizio fine], 'norm_factor', opt.Results.norm_factor, 'cache', 'motion_cache');
sm_freeze_full = sm_freeze_full .* ~is_below_threshold;

out = struct('bouts_proc', bouts_proc, 'sm_freeze_full', sm_freeze_full, ...
    'is_below_threshold', is_below_threshold, 'inizio', inizio, 'fine', fine);
end
