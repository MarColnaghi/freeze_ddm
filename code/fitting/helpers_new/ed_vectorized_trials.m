function result = ed_vectorized_trials(ts, theta, mu, fs, output_type)
    
    % Default to PDF output
    if nargin < 5
        output_type = 'pdf';
    end
    
    % Make sure mu is a matrix [n_trials × n_frames]
    if isvector(mu)
        mu = mu(:)';
    end
    [n_trials, n_frames] = size(mu);
    mu = mu .* (1/fs);

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
    
    % Time to frames (discrete)
    frames = floor(ts * fs);
    % frames = ts * fs;

    % Clip frames to valid range
    frames = max(1, min(n_frames, frames));
    
    % Compute all survival probabilities
    surv_probs = normcdf(theta - mu, 0, sqrt(1/fs));
    
    % Compute cumulative products (survival up to each frame)
    cum_surv = cumprod(surv_probs, 2);
    
    % Create linear indices for each trial's specific frame
    trial_indices = (1:n_trials)';
    
    % Compute CDF
    if strcmp(output_type, 'cdf') || strcmp(output_type, 'both')
        linear_idx = sub2ind([n_trials, n_frames], trial_indices, frames);
        survival_prob_at_t = cum_surv(linear_idx);
        cdf = survival_prob_at_t;
        cdf = cdf(:)'; % row vector
    end
    
    % Compute PDF
    if strcmp(output_type, 'pdf') || strcmp(output_type, 'both')
        % Get survival probability up to frame-1 for each trial
        survival_prob = ones(n_trials, 1);
        mask = frames > 1; % trials where frame > 1
        if any(mask)
            linear_idx = sub2ind([n_trials, n_frames], trial_indices(mask), frames(mask) - 1);
            survival_prob(mask) = cum_surv(linear_idx);
        end
        
        % Get mu value at the specific frame for each trial
        linear_idx_mu = sub2ind([n_trials, n_frames], trial_indices, frames);
        mu_at_frames = mu(linear_idx_mu);
        
        % Crossing probability at the specific frame for each trial
        cross_prob = 1 - normcdf(theta - mu_at_frames, 0, sqrt(1/fs));
        
        % Final probability [n_trials × 1]
        pdf = cross_prob .* survival_prob;
        pdf = pdf(:)'; % row vector
    end
    
    % Return based on output_type
    switch output_type
        case 'pdf'
            result = pdf;
        case 'cdf'
            result = cdf;
        case 'both'
            result.pdf = pdf;
            result.cdf = cdf;
        otherwise
            error('output_type must be ''pdf'', ''cdf'', or ''both''');
    end
end