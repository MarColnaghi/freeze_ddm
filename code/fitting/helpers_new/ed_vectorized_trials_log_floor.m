function result = ed_vectorized_trials_log_floor(ts, theta, mu, fs, output_type)
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

% Convert times to CONTINUOUS frames [n_trials × 1]
frames_continuous = ts * fs;

% Get integer part (j) and fractional part (alpha)
j = floor(frames_continuous);
alpha = frames_continuous - j;

% Clip j to valid range [1, n_frames-1] to ensure j+1 exists
j = max(1, min(n_frames - 1, j));

% Alpha should be in [0,1]
alpha = max(0, min(1, alpha));

% Compute survival probabilities [n_trials × n_frames]
surv_probs = normcdf(theta - mu, 0, sqrt(1/fs));

% Make sure that you don't get prob = 0 for log calculation
eps_prob = 1e-6;
surv_probs = max(eps_prob, min(1 - eps_prob, surv_probs));

% Compute log survival probabilities [n_trials × n_frames]
log_surv_probs = log(surv_probs);

% Compute cumulative sum in log space (survival up to each frame) [n_trials × n_frames]
log_cum_surv = cumsum(log_surv_probs, 2);

% Create linear indices for each trial
trial_indices = (1:n_trials)';

% Compute LOG CDF if needed
if strcmp(output_type, 'cdf') || strcmp(output_type, 'both')
    % Get log CDF at frame j and j+1
    linear_idx_j = sub2ind([n_trials, n_frames], trial_indices, j);
    linear_idx_j1 = sub2ind([n_trials, n_frames], trial_indices, j + 1);
    
    log_cdf_j = log_cum_surv(linear_idx_j);
    log_cdf_j1 = log_cum_surv(linear_idx_j1);
    
    % Interpolate in log space (actually need to be careful here)
    % For probabilities, better to interpolate in probability space then take log
    cdf_j = exp(log_cdf_j);
    cdf_j1 = exp(log_cdf_j1);
    
    % Linear interpolation: CDF = (1-alpha)*CDF_j + alpha*CDF_{j+1}
    cdf_interp = (1 - alpha) .* cdf_j + alpha .* cdf_j1;
    log_cdf = log(max(eps_prob, cdf_interp));
    
    % Replace -Inf with very negative number
    log_cdf(~isfinite(log_cdf)) = -200;
    log_cdf = log_cdf(:)'; % row vector
end

% Compute LOG PDF if needed
if strcmp(output_type, 'pdf') || strcmp(output_type, 'both')
    % Get mu values at frames j and j+1
    linear_idx_mu_j = sub2ind([n_trials, n_frames], trial_indices, j);
    linear_idx_mu_j1 = sub2ind([n_trials, n_frames], trial_indices, j + 1);
    
    mu_j = mu(linear_idx_mu_j);
    mu_j1 = mu(linear_idx_mu_j1);
    
    % Get LOG survival probability up to frame j-1 and j for interpolation
    log_survival_j = zeros(n_trials, 1);
    log_survival_j1 = zeros(n_trials, 1);
    
    mask_j = j > 1;
    if any(mask_j)
        linear_idx = sub2ind([n_trials, n_frames], trial_indices(mask_j), j(mask_j) - 1);
        log_survival_j(mask_j) = log_cum_surv(linear_idx);
    end
    
    % For j+1, survival up to frame j
    linear_idx = sub2ind([n_trials, n_frames], trial_indices, j);
    log_survival_j1 = log_cum_surv(linear_idx);
    
    % Crossing probability at frames j and j+1
    cross_prob_j = 1 - normcdf(theta - mu_j, 0, sqrt(1/fs));
    cross_prob_j1 = 1 - normcdf(theta - mu_j1, 0, sqrt(1/fs));
    
    cross_prob_j = max(eps_prob, min(1 - eps_prob, cross_prob_j));
    cross_prob_j1 = max(eps_prob, min(1 - eps_prob, cross_prob_j1));
    
    % Compute log PDF at j and j+1
    log_pdf_j = log(cross_prob_j) + log_survival_j;
    log_pdf_j1 = log(cross_prob_j1) + log_survival_j1;
    
    % Interpolate: LL = LL_j + alpha*(LL_{j+1} - LL_j)
    log_pdf = log_pdf_j + alpha .* (log_pdf_j1 - log_pdf_j);
    
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