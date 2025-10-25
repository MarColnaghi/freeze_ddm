function result = ed_vectorized_trials_log(ts, theta, mu, fs, output_type)

% Default to PDF output
if nargin < 5
    output_type = 'pdf';
end

% Make sure mu is a matrix [n_trials × n_frames]
if isvector(mu)
    mu = mu(:)';
end
[n_trials, n_frames] = size(mu);

% Make sure ts is column vector [n_trials × 1]
if isscalar(ts)
    ts = repmat(ts, n_trials, 1);
else
    ts = ts(:);
end

% Make sure theta is column vector [n_trials × 1]
if isscalar(theta)
    theta = repmat(theta, n_trials, 1);
else
    theta = theta(:);
end

% Check dimensions match
if length(ts) ~= n_trials || length(theta) ~= n_trials
    error('ts, theta, and number of rows in mu must match');
end

% Convert times to frames [n_trials × 1]
frames = round(ts * fs);

% Clip frames to valid range
frames = max(1, min(n_frames, frames));

% Compute survival probabilities [n_trials × n_frames]
surv_probs = normcdf(theta - mu, 0, sqrt(1/fs));

% Make sure that you don't get prob = 0 for log calculation
eps_prob = 1e-6; 
surv_probs = max(eps_prob, min(1 - eps_prob, surv_probs));

% Compute log survival probabilities [n_trials × n_frames]
log_surv_probs = log(surv_probs);

% Compute cumulative sum in log space (survival up to each frame) [n_trials × n_frames]
log_cum_surv = cumsum(log_surv_probs, 2);

% Create linear indices for each trial's specific frame
trial_indices = (1:n_trials)';

% Compute LOG CDF if needed
if strcmp(output_type, 'cdf') || strcmp(output_type, 'both')
    linear_idx = sub2ind([n_trials, n_frames], trial_indices, frames);
    log_survival_prob_at_t = log_cum_surv(linear_idx);

    % CDF = survival_prob_at_t (based on your original code)
    log_cdf = log_survival_prob_at_t;

    % Replace -Inf with very negative number
    log_cdf(~isfinite(log_cdf)) = -200;

    log_cdf = log_cdf(:)'; % row vector
end

% Compute LOG PDF if needed
if strcmp(output_type, 'pdf') || strcmp(output_type, 'both')
    % Get LOG survival probability up to frame-1 for each trial
    log_survival_prob = zeros(n_trials, 1);
    mask = frames > 1; % trials where frame > 1
    if any(mask)
        linear_idx = sub2ind([n_trials, n_frames], trial_indices(mask), frames(mask) - 1);
        log_survival_prob(mask) = log_cum_surv(linear_idx);
    end

    % Get mu value at the specific frame for each trial
    linear_idx_mu = sub2ind([n_trials, n_frames], trial_indices, frames);
    mu_at_frames = mu(linear_idx_mu);

    % Crossing probability at the specific frame for each trial
    cross_prob = 1 - normcdf(theta - mu_at_frames, 0, sqrt(1/fs));

    % Make sure that you don't get prob = 0 for log calculation
    cross_prob = max(eps_prob, min(1 - eps_prob, cross_prob));

    % LOG crossing probability
    log_cross_prob = log(cross_prob);

    % Final log probability [n_trials × 1]
    log_pdf = log_cross_prob + log_survival_prob;

    % Replace -Inf with very negative number
    log_pdf(~isfinite(log_pdf)) = -200;

    log_pdf = log_pdf(:)'; % row vector
end

% Return based on output_type
switch output_type
    case 'pdf'
        result = log_pdf;
    case 'cdf'
        result = log_cdf;
    case 'both'
        result.log_pdf = log_pdf;
        result.log_cdf = log_cdf;
    otherwise
        error('output_type must be ''pdf'', ''cdf'', or ''both''');
end
end