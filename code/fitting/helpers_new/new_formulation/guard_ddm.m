function y = guard_ddm(fun, t, mu, th, ndt)
% Returns 0 for entries where t <= ndt; calls `fun` otherwise
% fun  : function handle (e.g., pdf_raw.ddm)
% t    : vector of durations
% mu   : vector of drift rates
% th   : vector of thresholds
% ndt  : vector of non-decision times

    % Ensure all inputs are column vectors for consistency
    mu  = mu(:); 
    th  = th(:); 
    ndt = ndt(:);
    
    % Handle scalar t (like t0 or C) vs vector t
    if isscalar(t)
        t_vec = repmat(t, size(mu));
    else
        t_vec = t(:);
    end
    
    y = zeros(size(mu));
    
    % Only calculate for valid time points (t > ndt)
    mask = t_vec > ndt;
    if any(mask)
        % Call the raw DDM function only for the valid indices
        y(mask) = fun(t_vec(mask), mu(mask), th(mask), ndt(mask));
    end
end