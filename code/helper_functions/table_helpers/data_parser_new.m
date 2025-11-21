
function [temp] = data_parser_new(bouts, varargin)

%% Select immobility bouts and Rescale the Variables for Model Fitting

opt = inputParser;
addParameter(opt, 'type', 'immobility');
addParameter(opt, 'period', 'loom');
addParameter(opt, 'window', '');
addParameter(opt, 'nloom', []);
addParameter(opt, 'sloom', []);
addParameter(opt, 'min_dur', 0);
parse(opt, varargin{:});

% Select immobility bouts
if strcmp(opt.Results.type, 'immobility')
    temp = bouts(bouts.type == 1,:);
elseif strcmp(opt.Results.type, 'mobility')
    temp = bouts(bouts.type == 0,:);
else
    temp = bouts;
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
    temp = temp(temp.nloom == [opt.Results.nloom], :);
end

if ~isempty(opt.Results.sloom)
    temp = temp(temp.sloom == [opt.Results.sloom], :);
end

% Additional filters
temp = temp(temp.durations >= opt.Results.min_dur, :);

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
