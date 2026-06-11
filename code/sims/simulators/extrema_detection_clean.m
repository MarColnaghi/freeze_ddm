function [rt_frames, traj, t] = extrema_detection_clean(signal, theta, sigma, dt, seed)
% EXTREMA_DETECTION_CLEAN  First-passage of an instantaneous / extrema detector.
%
%   [rt_frames, traj, t] = extrema_detection_clean(signal, theta, sigma, dt, seed)
%
%   Un-freeze at the FIRST frame whose momentary (noisy) signal exceeds theta:
%       m(k) = signal(k) + sigma*sqrt(dt)*randn
%       decision: first k with m(k) >= theta
%   There is NO accumulation — only the current frame's value matters. The
%   recovered hazard kernel should therefore be a DELTA at tau = 0.
%
%   signal : per-frame instantaneous drive LEVEL (vector), e.g. mu0 + beta*sm(t).
%            (Note: here drive is a level, not a rate — different units from the
%             accumulator's drift; tune theta/sigma accordingly.)
%   theta  : decision bound.   sigma : diffusion SD.   dt : seconds/frame (1/fps).
%
%   Outputs:
%       rt_frames : 1-based frame index of first crossing, or NaN if none.
%       traj      : momentary value m(k) at every frame (column, length numel(signal));
%                   NaN for frames after the crossing.
%       t         : time vector, t(k) = k*dt, so rt_frames*dt = t(rt_frames).
%
%   The noise term (sigma*sqrt(dt)) and crossing rule (>= theta) are now SHARED
%   with sim_leaky_accumulator / the C++ mex, so the momentary detector and the
%   accumulators see the same per-frame noise and their trajectories are
%   directly comparable under a common seed (both draw one randn per frame, in
%   frame order). Note this makes the effective noise dt-dependent — to keep a
%   dt-robust momentary rule instead, scale sigma by 1/sqrt(dt) at the call site.
%
%   This replaces the old extrema_detection_new.m, which (a) had a length bug
%   ([z; n_t-vector] -> n_t+1 mismatched with t), (b) crossed on the dt-scaled
%   per-step INCREMENT (not a dt-robust instantaneous rule), and (c) carried a
%   wrong (copy-pasted DDM) docstring.

    if nargin >= 5 && ~isempty(seed), rng(seed); end
    n      = numel(signal);
    nscale = sigma*sqrt(dt);
    t      = (1:n)' * dt;
    traj   = nan(n, 1);
    rt_frames = NaN;
    for k = 1:n
        m = signal(k) + nscale*randn;
        traj(k) = m;
        if m >= theta
            rt_frames = k;
            return
        end
    end
end
