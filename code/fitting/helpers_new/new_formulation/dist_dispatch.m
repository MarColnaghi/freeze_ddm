function [f, F] = dist_dispatch(type, out, prefix, ndt, t)
    % This handles the vector math for each process type.
    % 'out' is the evaluated parameter struct from evaluate_model_fast
    
    % Ensure t is treated as a vector if it's a scalar constant (t0/C)
    if isscalar(t)
        t = repmat(t, size(ndt)); 
    end
    
    switch type
        case 'ddm'
            [pdf_raw, cdf_raw] = pdf_cdf({'ddm'});
            % Use your existing guard_ddm to handle t < ndt
            f = guard_ddm(pdf_raw.ddm, t, out.([prefix '_mu']), out.([prefix '_theta']), ndt);
            F = guard_ddm(cdf_raw.ddm, t, out.([prefix '_mu']), out.([prefix '_theta']), ndt);
            
        case 'exp'
            [pdf_raw, cdf_raw] = pdf_cdf({'exp'});
            % Standard Exponential: f(t) = lambda * exp(-lambda * (t-ndt))
            % Assuming 'mu' was mapped to the rate parameter
            lambda = out.([prefix '_mu']);
            dt = max(0, t - ndt);
            f = lambda .* exp(-lambda .* dt);
            F = 1 - exp(-lambda .* dt);
            
        case 'ed'
            % Add your specific Exponential Decay logic here
            f = zeros(size(t)); F = zeros(size(t)); 
    end
end

