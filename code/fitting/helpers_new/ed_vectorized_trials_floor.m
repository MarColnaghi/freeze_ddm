function result = ed_vectorized_trials_floor(ts, theta, mu, fs, output_type)
% Non-log space version with interpolation
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

% Clamp probabilities to avoid numerical issues
eps_prob = 1e-10;
surv_probs = max(eps_prob, min(1 - eps_prob, surv_probs));

% Compute cumulative product (survival up to each frame) [n_trials × n_frames]
cum_surv = cumprod(surv_probs, 2);

% Create linear indices for each trial
trial_indices = (1:n_trials)';

% Compute CDF if needed
if strcmp(output_type, 'cdf') || strcmp(output_type, 'both')
    % Get survival at frame j and j+1
    linear_idx_j = sub2ind([n_trials, n_frames], trial_indices, j);
    linear_idx_j1 = sub2ind([n_trials, n_frames], trial_indices, j + 1);
    
    surv_j = cum_surv(linear_idx_j);
    surv_j1 = cum_surv(linear_idx_j1);
    
    % Interpolate survival in GEOMETRIC space
    % S(t) = S_j^(1-alpha) * S_{j+1}^alpha
    surv_interp = surv_j.^(1 - alpha) .* surv_j1.^alpha;
    surv_interp = max(eps_prob, min(1 - eps_prob, surv_interp));
    
    % IMPORTANT: Based on your working code, "cdf" actually means survival!
    cdf_result = surv_interp(:)'; % row vector (this is actually survival)
end

% Compute PDF if needed
if strcmp(output_type, 'pdf') || strcmp(output_type, 'both')
    % Interpolate mu values at the fractional frame position
    linear_idx_mu_j = sub2ind([n_trials, n_frames], trial_indices, j);
    linear_idx_mu_j1 = sub2ind([n_trials, n_frames], trial_indices, j + 1);
    
    mu_j = mu(linear_idx_mu_j);
    mu_j1 = mu(linear_idx_mu_j1);
    
    % Linear interpolation of mu at fractional frame
    mu_interp = (1 - alpha) .* mu_j + alpha .* mu_j1;
    
    % Get survival probability up to BEFORE the fractional frame
    % Use geometric interpolation of survival
    survival_before_j = ones(n_trials, 1);
    survival_before_j1 = ones(n_trials, 1);
    
    mask_j = j > 1;
    if any(mask_j)
        linear_idx = sub2ind([n_trials, n_frames], trial_indices(mask_j), j(mask_j) - 1);
        survival_before_j(mask_j) = cum_surv(linear_idx);
    end
    
    linear_idx = sub2ind([n_trials, n_frames], trial_indices, j);
    survival_before_j1 = cum_surv(linear_idx);
    
    % Interpolate survival up to before crossing
    survival_before_interp = survival_before_j.^(1 - alpha) .* survival_before_j1.^alpha;
    
    % Crossing probability at the interpolated mu
    cross_prob_interp = 1 - normcdf(theta - mu_interp, 0, sqrt(1/fs));
    cross_prob_interp = max(eps_prob, min(1 - eps_prob, cross_prob_interp));
    
    % PDF = P(cross at interpolated t) * P(survived up to interpolated t)
    pdf_interp = cross_prob_interp .* survival_before_interp;
    pdf_interp = max(eps_prob, pdf_interp);
    
    pdf_result = pdf_interp(:)'; % row vector
end

% Return based on output_type
switch output_type
    case 'pdf'
        result = pdf_result;
    case 'cdf'
        result = cdf_result;
    case 'both'
        result.pdf = pdf_result;
        result.cdf = cdf_result;
    otherwise
        error('output_type must be ''pdf'', ''cdf'', or ''both''');
end
end