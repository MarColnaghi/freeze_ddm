function prob = ed_at_time(t_sec, theta, mu, fs)
    % Convert seconds to frame index
    frame = round(t_sec * fs);
    
    % Check bounds
    if frame < 1 || frame > length(mu)
        prob = 0;
        return;
    end
    
    % Probability of crossing at this frame
    cross_prob = 1 - normcdf(theta - mu(frame));
    
    % Probability of NOT crossing at all previous frames
    if frame > 1
        survival_prob = prod(normcdf(theta - mu(1:frame-1)));
    else
        survival_prob = 1;
    end
    
    prob = cross_prob * survival_prob;
end