function [rt, traj, t, det_cum, sto_cum] = extrema_detection_new(varargin)

%
%   [rt, traj, t] = extrema_detection_new('mu', val, 'mu_t', vec, ...)
%   simulates a DDM with optional time-varying drift.
%
%   Parameters:
%       'mu'       - Constant drift rate (ignored if 'mu_t' provided)
%       'mu_t'     - Time-varying drift rate (vector)
%       'theta'    - Decision threshold
%       'z'        - Starting point
%       'sigma'    - Noise standard deviation
%       'dt'       - Time step
%       'T'        - Max simulation time
%       'ndt'      - Non-decision time
%       'truncate' - (true/false) restrict RT to [0.5, 10.7]
%
%   Outputs:
%       rt   - Reaction time
%       traj - Trajectory (NaNs after decision)
%       t    - Time vector

    p = inputParser;
    addParameter(p, 'mu',       [], @(x) validateattributes(x, {'numeric'}, {'scalar'}));
    addParameter(p, 'mu_t',     [], @(x) validateattributes(x, {'numeric'}, {'vector'}));
    addParameter(p, 'theta',    [], @(x) validateattributes(x, {'numeric'}, {'scalar'}));
    addParameter(p, 'z',        0,  @(x) validateattributes(x, {'numeric'}, {'scalar'}));
    addParameter(p, 'sigma',    1,  @(x) validateattributes(x, {'numeric'}, {'scalar', 'nonnegative'}));
    addParameter(p, 'dt',       0.01, @(x) validateattributes(x, {'numeric'}, {'scalar', 'positive'}));
    addParameter(p, 'T',        10,   @(x) validateattributes(x, {'numeric'}, {'scalar', 'positive'}));
    addParameter(p, 'ndt',      0.3,  @(x) validateattributes(x, {'numeric'}, {'scalar', 'nonnegative'}));
    addParameter(p, 'truncate', false, @(x) islogical(x));
    addParameter(p, 'seed', []);

    parse(p, varargin{:});
    params = p.Results;

    if ~isempty(params.seed)
        rng(params.seed);
    end

    % Time vector
    t = 0:params.dt:params.T;
    n_t = numel(t);

    % Use time-varying drift if provided
    use_dynamic_mu = ~isempty(params.mu_t);
    if use_dynamic_mu
        if length(params.mu_t) ~= n_t - 1
            error('mu_t must match the number of time steps.');
        end
        mu_t = params.mu_t(:);  % Ensure column vector
    else
        mu_t = params.mu * ones(n_t - 1, 1);
    end

    while true
         dW = randn(n_t - 1, 1);
        det_part = mu_t * params.dt ;
        sto_part = params.sigma * sqrt(params.dt) .* dW;
        dx = [params.z; det_part + sto_part];
        det_cum = [params.z; cumsum(det_part)];
        sto_cum = [0; cumsum(sto_part)];
        traj = dx;
        rt = NaN;

        % Check threshold crossing
        cross_idx = find(traj >= params.theta, 1);
        if ~isempty(cross_idx)
            rt = t(cross_idx) + params.ndt;
            traj((cross_idx+1):end) = NaN;

        end

        if ~params.truncate || (~isnan(rt) && rt > 0.5 && rt <= 10.7)
            break;
        end
    end
end
