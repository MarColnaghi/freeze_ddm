function [rt_frames, traj, t] = sim_leaky_accumulator(drift, theta, sigma, lambda, dt, seed, x0)
% SIM_LEAKY_ACCUMULATOR  First-passage time of a leaky (or perfect) accumulator.
%
%   [rt_frames, traj, t] = sim_leaky_accumulator(drift, theta, sigma, lambda, dt, seed, x0)
%
%   Euler–Maruyama integration of
%       x(k+1) = x(k)*(1 - lambda*dt) + drift(k)*dt + sigma*sqrt(dt)*randn
%   starting from x = x0 (default 0), returning the 1-based frame index of the
%   first time x >= theta, or NaN if it never crosses within numel(drift) steps.
%
%   x0 : initial accumulator position (optional, default 0). ABSOLUTE units,
%        same scale as theta (matches params[2] of the sim_ddm_seeded mex). Pass
%        e.g. initial_value_frac*theta for a bayes_fpe-style start point.
%
%   This matches the C++ mex (sim_ddm_seeded / sim_ddm_ultra) exactly:
%     lambda = 0      -> PERFECT accumulator (DDM; equal weighting of past drift)
%     lambda = 1/tau  -> LEAKY integrator with time-constant tau (recency-weighted)
%
%   drift : per-frame drift RATE (vector), e.g. mu0 + beta*sm(t).
%   theta : decision bound.   sigma : diffusion SD.   dt : seconds/frame (1/fps).
%
%   Outputs:
%       rt_frames : 1-based frame index of first crossing, or NaN if none.
%       traj      : accumulator state x at every frame (column, length numel(drift));
%                   NaN for frames after the crossing. traj(k) is the state AFTER
%                   integrating frame k.
%       t         : time vector, t(k) = k*dt, so rt_frames*dt = t(rt_frames).
%
%   The leaky integrator's first-passage hazard depends on a leaky integral of
%   the drift with exponential memory exp(-lambda*tau) -> the recovered kernel
%   should be exponential with timescale tau. lambda=0 -> flat kernel.
%
%   Noise (sigma*sqrt(dt)) and crossing rule (x >= theta) are shared with
%   extrema_detection_clean so trajectories are directly comparable under a
%   common seed (both draw one randn per frame, in frame order).

    if nargin >= 6 && ~isempty(seed), rng(seed); end
    if nargin < 7 || isempty(x0), x0 = 0; end
    n      = numel(drift);
    decay  = 1 - lambda*dt;
    nscale = sigma*sqrt(dt);   % seam: a per-frame sigma vector (signal-dependent
                               % noise) would make this an n-vector indexed by k.
    t      = (1:n)' * dt;
    traj   = nan(n, 1);
    x      = x0;
    rt_frames = NaN;
    for k = 1:n
        x = x*decay + drift(k)*dt + nscale*randn;
        traj(k) = x;
        if x >= theta
            rt_frames = k;
            return
        end
    end
end
