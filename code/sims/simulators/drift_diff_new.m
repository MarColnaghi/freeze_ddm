function [rt, traj, t, det_cum, sto_cum] = drift_diff_new(varargin)
    p = inputParser;
    addParameter(p, 'mu',       0,    @(x) validateattributes(x, {'numeric'}, {'scalar'}));
    addParameter(p, 'mu_t',     [],   @(x) validateattributes(x, {'numeric'}, {'vector'}));
    addParameter(p, 'theta',    1,    @(x) validateattributes(x, {'numeric'}, {'scalar'}));
    addParameter(p, 'z',        0,    @(x) validateattributes(x, {'numeric'}, {'scalar'}));
    addParameter(p, 'sigma',    1,    @(x) validateattributes(x, {'numeric'}, {'scalar', 'nonnegative'}));
    addParameter(p, 'dt',       0.01, @(x) validateattributes(x, {'numeric'}, {'scalar', 'positive'}));
    addParameter(p, 'T',        10.5, @(x) validateattributes(x, {'numeric'}, {'scalar', 'positive'}));
    addParameter(p, 'ndt',      0,    @(x) validateattributes(x, {'numeric'}, {'scalar', 'nonnegative'}));
    addParameter(p, 'truncate', false,@(x) islogical(x));
    addParameter(p, 'seed',     []);
    
    parse(p, varargin{:});
    params = p.Results;
    if ~isempty(params.seed), rng(params.seed); end
    
    t = 0:params.dt:params.T;
    n_t = numel(t);
    n_steps = n_t - 1;

    if ~isempty(params.mu_t)
        if length(params.mu_t) ~= n_steps
            error('mu_t must have length %d (n_t - 1).', n_steps);
        end
        mu_t = params.mu_t(:);
    else
        mu_t = params.mu * ones(n_steps, 1);
    end

    while true
        dW = randn(n_steps, 1);
        
        % 1. Calculate the raw displacements (deltas)
        det_part_raw = mu_t * params.dt;
        sto_part_raw = params.sigma * sqrt(params.dt) .* dW;
        
        % 2. INTEGRATION: traj = z + sum of all previous steps
        % We use [0; cumsum(...)] so that the first point is ALWAYS exactly 'z'
        det_cum = [params.z; params.z + cumsum(det_part_raw)];
        sto_cum = [0; cumsum(sto_part_raw)];
        traj = det_cum + sto_cum;
        
        rt = NaN;
        cross_idx = find(traj >= params.theta, 1);
        
        if ~isempty(cross_idx)
            rt = t(cross_idx) + params.ndt;
            % Apply the NaNs for visualization/output
            traj(cross_idx+1:end) = NaN;
            det_cum(cross_idx+1:end) = NaN;
            sto_cum(cross_idx+1:end) = NaN;
        end
        
        if ~params.truncate || (~isnan(rt) && rt > 0.5 && rt <= 10.7)
            break;
        end
    end
end