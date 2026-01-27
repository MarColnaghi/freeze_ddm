function result = ed_vectorized_trials(ts, theta, mu, fs, output_type)
% ED_VECTORIZED_TRIALS
% Discrete-time first-passage PDF / CDF for an evidence diffusion process
% with time-varying drift.
%
% ts          : [n_trials x 1] crossing times (s)
% theta       : [n_trials x 1] decision thresholds
% mu          : [n_trials x n_frames] drift per frame
% fs          : sampling frequency (Hz)
% output_type : 'pdf', 'cdf', 'both', 'logpdf', 'logcdf'
%
% Notes:
% - Continuous times are discretized using ceil(ts * fs)
% - Computations are done in log-space for numerical stability

    % -------------------------------
    % Defaults
    % -------------------------------
    if nargin < 5
        output_type = 'pdf';
    end

    % -------------------------------
    % Shape handling
    % -------------------------------
    if isvector(mu)
        mu = mu(:)';
    end
    [n_trials, n_frames] = size(mu);

    % Scale drift to per-frame units
    mu = mu ./ fs;

    % ts
    if isscalar(ts)
        ts = repmat(ts, n_trials, 1);
    else
        ts = ts(:);
    end

    % theta
    if isscalar(theta)
        theta = repmat(theta, n_trials, 1);
    else
        theta = theta(:);
    end

    if length(ts) ~= n_trials || length(theta) ~= n_trials
        error('ts, theta, and number of rows in mu must match');
    end

    % -------------------------------
    % Time → frame discretization
    % -------------------------------
    % Trials at or before t = 0
    too_early = ts <= 0;

    % Discretize time (first frame >= ts)
    frames = ceil(ts * fs);

    % Clip to valid range
    frames = max(1, min(n_frames, frames));

    % -------------------------------
    % Survival probabilities (log-space)
    % -------------------------------
    sigma = sqrt(1/fs);

    % Per-frame survival prob
    surv_probs = normcdf(theta - mu, 0, sigma);

    % Avoid log(0)
    surv_probs = max(surv_probs, realmin);

    % Log survival
    log_surv = log(surv_probs);

    % Cumulative log survival
    cum_log_surv = cumsum(log_surv, 2);

    trial_idx = (1:n_trials)';

    % -------------------------------
    % CDF
    % -------------------------------
    if any(strcmp(output_type, {'cdf','both','logcdf'}))
        lin_idx = sub2ind([n_trials, n_frames], trial_idx, frames);
        log_survival_at_t = cum_log_surv(lin_idx);

        logcdf = log_survival_at_t;
        logcdf(too_early) = 0;   % log(1)

        cdf = exp(logcdf);
    end

    % -------------------------------
    % PDF
    % -------------------------------
    if any(strcmp(output_type, {'pdf','both','logpdf'}))
        % Survival up to frame-1
        log_survival_prev = zeros(n_trials,1); % log(1)

        mask = frames > 1;
        if any(mask)
            lin_idx_prev = sub2ind([n_trials, n_frames], ...
                                   trial_idx(mask), frames(mask)-1);
            log_survival_prev(mask) = cum_log_surv(lin_idx_prev);
        end

        % Survival at frame (reuse!)
        lin_idx_mu = sub2ind([n_trials, n_frames], trial_idx, frames);
        surv_at_frame = surv_probs(lin_idx_mu);

        % Crossing probability
        cross_prob = max(1 - surv_at_frame, realmin);

        logpdf = log(cross_prob) + log_survival_prev;
        logpdf(too_early) = -Inf;

        pdf = exp(logpdf);
    end

    % -------------------------------
    % Output
    % -------------------------------
    switch output_type
        case 'pdf'
            result = pdf';
        case 'cdf'
            result = cdf';
        case 'logpdf'
            result = logpdf';
        case 'logcdf'
            result = logcdf';
        case 'both'
            result.pdf    = pdf';
            result.cdf    = cdf';
            result.logpdf = logpdf';
            result.logcdf = logcdf';
        otherwise
            error('Invalid output_type');
    end
end
