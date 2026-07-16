function [rt, info] = simulate_freeze_ddm(signal, p, N)
% SIMULATE_FREEZE_DDM  Forward-simulate first-passage times with the full
% bayes_fpe-style parameter set, reusing the existing freeze_ddm accumulators.
%
%   [rt, info] = simulate_freeze_ddm(signal, p, N)
%
%   Builds the effective per-frame drift with SIGNAL_TO_DRIFT (ReLU, filter,
%   clip, sensory delay, gain/bias -- none of which touch the dynamics), then
%   runs N Monte-Carlo trials through a leaky accumulator. The accumulator-level
%   parameters (leak, bound, initial position) and the observation-level shifts
%   (delayed start, non-decision time, contaminant lapse) are applied here.
%
%   signal : raw per-frame signal for ONE trial (vector, t=0 at freeze onset).
%   N      : number of Monte-Carlo trials.
%   p      : parameter struct. Fields (all optional unless noted):
%     dt                 seconds/frame (REQUIRED).
%     theta              decision bound (default 1). For loom/focal bound
%                        structure, compute theta per trial in the caller and
%                        pass the scalar here.
%     sigma              diffusion SD (default 1). NB the mex takes sigma^2;
%                        this wrapper squares it for you.
%     leak               leak rate lambda, 1/s (default 0 = perfect integrator).
%     initial_value_frac start point as a fraction of theta (default 0);
%                        x0 = initial_value_frac*theta (absolute units).
%     delayed_start      accumulation start delay, s (default 0): the first
%                        round(delayed_start/dt) frames are skipped and the
%                        value is added back to the returned RT.
%     ndt                non-decision time, s (default 0), added to the RT.
%     contaminant_prob   lapse probability (default 0): with this prob the RT is
%                        replaced by Uniform(0, t_max).
%     backend            'mex' | 'matlab' | 'auto' (default 'auto' -> 'mex').
%     seed               base RNG seed (optional; reproducible if set).
%     ...plus every SIGNAL_TO_DRIFT field (gain, bias, relu_*, filter_tau,
%        clip_threshold, sensory_delay, time_kernel_*).
%
%   rt   : N x 1 first-passage times in seconds (NaN = no crossing / censored).
%   info : struct with the resolved drift, theta, x0, backend, and t_max.
%
%   NOTE (known divergence from the fitted FPE model): the freeze_ddm
%   accumulator has only an absorbing UPPER bound -- there is no reflecting
%   lower boundary like bayes_fpe's x_min. Under strong leak the low end of the
%   trajectory can therefore diverge from the FPE model. Signal-dependent
%   diffusion noise (signal_noise_scale) is intentionally NOT supported here.

    if nargin < 3 || isempty(N), N = 1; end
    dt = must(p, 'dt');

    theta   = getf(p, 'theta', 1);
    sigma   = getf(p, 'sigma', 1);
    lambda  = getf(p, 'leak', 0);
    ivf     = getf(p, 'initial_value_frac', 0);
    dstart  = getf(p, 'delayed_start', 0);
    ndt     = getf(p, 'ndt', 0);
    pcont   = getf(p, 'contaminant_prob', 0);
    backend = lower(string(getf(p, 'backend', 'auto')));
    seed    = getf(p, 'seed', []);

    % ---- form the effective drift (all pre-accumulation preprocessing) ------
    drift = signal_to_drift(signal, p);
    drift = drift(:);
    t_max = numel(drift) * dt;                 % full pre-slice grid, for lapses

    % ---- start-timing: skip the pre-start frames, start at x0 = frac*theta --
    n_skip = round(dstart / dt);
    n_skip = max(0, min(n_skip, numel(drift) - 1));
    drift_sim = drift(n_skip + 1 : end);
    x0 = ivf * theta;

    if backend == "auto", backend = "mex"; end

    % ---- run N trials -------------------------------------------------------
    switch backend
        case "mex"
            % params = [dt, sigma^2, x0, theta, lambda]; one call -> N trials.
            mex_seed = seed;
            if isempty(mex_seed), mex_seed = randi(2^31 - 1); end
            fpt = sim_ddm_seeded(drift_sim, [dt, sigma^2, x0, theta, lambda], N, mex_seed);
            fpt = fpt(:);
        case "matlab"
            fpt = nan(N, 1);
            for i = 1:N
                s_i = [];
                if ~isempty(seed), s_i = seed + i; end     % [] -> independent draws
                k = sim_leaky_accumulator(drift_sim, theta, sigma, lambda, dt, s_i, x0);
                if ~isnan(k), fpt(i) = k * dt; end
            end
        otherwise
            error('simulate_freeze_ddm: unknown backend "%s".', backend);
    end

    % ---- observation-level: add back the start delay + ndt ------------------
    rt = fpt + dstart + ndt;                   % NaN stays NaN (censored)

    % ---- contaminant lapses: Uniform(0, t_max) ------------------------------
    if pcont > 0
        is_lapse = rand(N, 1) < pcont;
        rt(is_lapse) = rand(sum(is_lapse), 1) * t_max;
    end

    if nargout > 1
        info = struct('drift', drift, 'drift_sim', drift_sim, 'theta', theta, ...
                      'x0', x0, 'n_skip', n_skip, 'backend', char(backend), ...
                      't_max', t_max);
    end
end

% =========================================================================
function v = getf(p, name, default)
    if isstruct(p) && isfield(p, name) && ~isempty(p.(name))
        v = p.(name);
    else
        v = default;
    end
end

function v = must(p, name)
    assert(isstruct(p) && isfield(p, name) && ~isempty(p.(name)), ...
           'simulate_freeze_ddm: p.%s is required.', name);
    v = p.(name);
end
