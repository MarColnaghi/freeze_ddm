function [early_sm, late_sm, com, slopes, slopes_smoothed, com_norm, trend_score, diff_start_end] = partition_variable_freezes(mcell, varargin)

% INPUT:
%   mcell - 1 x N cell array, each cell contains a vector of motion values (length varies)
%   opt - 

% OUTPUT:
%   early_sm - indices of bouts with early CoM and decreasing motion
%   late_sm  - indices of bouts with late CoM and increasing motion
%   com     - center of mass values for each bout
%   slopes        - linear trend (slope) for each bout
%   com_norm       - normalized center of mass [0, 1]



    opt = inputParser;
    addParameter(opt, 'thr_com', 0.35);
    addParameter(opt, 'thr_slope', 0);
    addParameter(opt, 'thr_trend', 0);
    addParameter(opt, 'mov_mean', 12);
    addParameter(opt, 'thr_es', 0);

    parse(opt, varargin{:});

    N = numel(mcell);
    com = nan(1, N);
    com_norm = nan(1, N);
    slopes = nan(1, N);
    slopes_smoothed = nan(1, N);
    trend_score = nan(1, N);
    diff_start_end = nan(1, N);

    for i = 1:N
        m = mcell{i}(1:end);
        T = length(m);
        t = (1:T)';
        smoothed = movmean(m(:), opt.Results.mov_mean);

        if sum(m) == 0
            continue; 
        end

        % Center of Mass
        com(i) = sum(t .* smoothed(:)) / sum(smoothed);
        com_norm(i) = (com(i) - 1) / (T - 1);

        % End - Beginning difference
        w = max(1, floor(length(smoothed) * 0.2));  % Use 20% of the time series, or at least 1 point
        diff_start_end(i) = (mean(smoothed(end - w + 1:end)) - mean(smoothed(1:w)));% ./ (length(smoothed) * 0.2);

        % Trend (linear slope)
        p = polyfit(t, m(:), 1);
        slopes(i) = p(1);

        p = polyfit(1:length(smoothed), smoothed, 1);
        slopes_smoothed(i) = p(1);

        steps = diff(smoothed(:));
        trend_score(i) = sum(steps > 0) / length(steps);
    end

    % Define partitions
    early_sm = find(com_norm < opt.Results.thr_com & slopes <= -opt.Results.thr_slope & trend_score < opt.Results.thr_trend & diff_start_end < -opt.Results.thr_es);
    late_sm  = find(com_norm > (1 - opt.Results.thr_com) & slopes >= opt.Results.thr_slope & trend_score > (1 - opt.Results.thr_trend) & diff_start_end > opt.Results.thr_es);
end
