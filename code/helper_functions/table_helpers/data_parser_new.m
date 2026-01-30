
function [temp] = data_parser_new(bouts, varargin)

%% Select immobility bouts and Rescale the Variables for Model Fitting

opt = inputParser;
addParameter(opt, 'type', 'immobility');
addParameter(opt, 'period', '');
addParameter(opt, 'window', '');
addParameter(opt, 'nloom', []);
addParameter(opt, 'sloom', []);
addParameter(opt, 'min_dur', 0);
addParameter(opt, 'exclude_flies', false);

parse(opt, varargin{:});

% Select immobility bouts
if strcmp(opt.Results.type, 'immobility')
    temp = bouts(bouts.type == 1,:);
elseif strcmp(opt.Results.type, 'mobility')
    temp = bouts(bouts.type == 0,:);
else
    temp = bouts;
end

% Additional filters
temp = temp(temp.durations >= opt.Results.min_dur, :);

% Calculate the number of freezes per each fly.
[G, ~] = findgroups(temp.fly);
cell_array = splitapply(@(x) {(1:numel(x))'}, temp.fly, G);
temp.n_generated_freezes = vertcat(cell_array{:});

% Time since last time freeze
% Define a function that calculates Diff(Start_i - End_i-1)
% The first element is always NaN
calc_interval = @(o, e) {[o(1); o(2:end) - e(1:end-1)]};
time_since_last = splitapply(calc_interval, temp.onsets, temp.ends, G);
temp.time_since_last = vertcat(time_since_last{:});

% Calculate cumulative time spent freezing
cum_freeze_time = splitapply(@(x) {cumsum([0; x(1:end - 1)])}, temp.durations, G);
temp.cum_freeze_time = vertcat(cum_freeze_time{:});

% 1. Duration of the previous freeze
% Logic: Shift durations by 1. First element is NaN.
calc_prev_dur = @(x) {[NaN; x(1:end-1)]};
prev_dur = splitapply(calc_prev_dur, temp.durations, G);
temp.prev_freeze_dur = vertcat(prev_dur{:});

% 2. Average duration of freezes up until that freeze (Running Mean of History)
% Logic: Cumulative Sum of previous durations / Count of previous durations
% The first element is NaN because there is no history.
calc_avg_hist = @(x) {[NaN; cumsum(x(1:end-1)) ./ (1:numel(x)-1)']};
avg_history = splitapply(calc_avg_hist, temp.durations, G);
temp.avg_history_dur = vertcat(avg_history{:});

if opt.Results.exclude_flies
    % Only baseline immobility bouts
    bsl = bouts(bouts.type == 1 & bouts.period == 0, :);
    [G, fly_ids] = findgroups(bsl.fly);
    total_bsl_freeze = splitapply(@sum, bsl.durations, G);  % in frames
    thr = 150 * 60;  % e.g. 10 seconds in frames
    bad_flies = fly_ids(total_bsl_freeze > thr);
    temp = temp(~ismember(temp.fly, bad_flies), :);
    fprintf('Excluded %d flies due to high baseline freezing\n', bad_flies);
end

% Select period: bsl or loom
if strcmp(opt.Results.period, 'bsl')
    temp = temp(temp.period == 0, :);
    temp = temp(temp.nloom < 21, :);
    temp = temp(temp.frozen_start == 0, :);

elseif strcmp(opt.Results.period, 'loom')
    temp = temp(temp.period == 1, :);
    temp = temp(temp.frozen_start == 0, :);
end

% Select Window: le or not?
if strcmp(opt.Results.window, 'le')
    temp = temp(temp.le == 1, :);
end

% Select Loom Number
if ~isempty(opt.Results.nloom)
    temp = temp(ismember(temp.nloom, [opt.Results.nloom]), :);
end

if ~isempty(opt.Results.sloom)
    temp = temp(temp.sloom == [opt.Results.sloom], :);
end

% Rescale all the Variables
temp.durations_s = temp.durations/60;
temp.latency_s = temp.onsets_loomaligned/60;
temp.avg_sm_freeze_norm = temp.avg_sm/10;
temp.avg_fs_freeze_norm = temp.avg_fs/10;
temp.avg_fs_1s_norm = temp.avg_fs_1s/10;
temp.sloom_norm = temp.sloom/25;
temp.nloom_norm = temp.nloom/10;
temp.onsets_loomaligned_norm = temp.onsets_loomaligned/20;

% Set Thresholds on Freeze Durations, Social Motion and Speed of the Fly
motion_max = 2; % outliers
speed_max = 2; % outliers

% Threshold the Dataset

temp = temp(temp.avg_sm_freeze_norm <= motion_max & temp.avg_fs_1s_norm <= speed_max, :);


end
