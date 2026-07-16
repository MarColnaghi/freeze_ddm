function [drift, full] = signal_to_drift(signal, p)
% SIGNAL_TO_DRIFT  Turn a raw per-frame signal into an effective drift vector.
%
%   [drift, full] = signal_to_drift(signal, p)
%
%   Ports the bayes_fpe_fork signal-preprocessing front-end so the freeze_ddm
%   Monte-Carlo simulators can use the same knobs the fitted (Fokker-Planck)
%   model uses. NONE of these operations touch the accumulator dynamics: they
%   only reshape the per-frame drive BEFORE it is integrated. The accumulator
%   parameters (leak, bound, initial position) are handled in the caller
%   (see SIMULATE_FREEZE_DDM), not here.
%
%   signal : per-frame signal. A vector (1xT or Tx1) for one trial, or an
%            [nTrials x T] matrix (operated on row-wise). Output matches the
%            input orientation/shape.
%   p      : struct of parameters. Every field is optional; an unset field is a
%            no-op, so with defaults this reduces to  drift = bias + gain*signal
%            (i.e. the current compare_simulators.m:43 behaviour with bias=0).
%
%   Pipeline order (mirrors tv_core.py _apply_shared_signal_transforms, then the
%   time kernel; see bayes_fpe_fork/model_helpers/tv_core.py):
%     1. sensory delay        (time-shift)      p.sensory_delay, p.sensory_delay_mode
%     2. causal low-pass       (filter)          p.filter_tau
%     3. saturating clip                         p.clip_threshold
%     4. ReLU (thr/offset/pow)                   p.relu_threshold/relu_offset/relu_power
%     5. sigmoid time kernel                     p.time_kernel_midpoint/time_kernel_slope
%     6. drift = bias + gain .* processed        p.bias, p.gain
%
%   p.dt (seconds/frame) is required whenever a time-based step is used
%   (sensory delay, filter, or time kernel).
%
%   PRE-ONSET STREAM (p.n_pre, default 0)
%   -------------------------------------
%   A sensory_delay L means that at freeze-relative time tau the accumulator uses
%   evidence from raw[tau - L]. For tau < L that index is NEGATIVE, i.e. BEFORE
%   freeze onset. bayes_fpe fills that window with ZEROS -- not a physiological
%   claim, just a consequence of its signal array starting at freeze onset
%   (delay_tv_signals_with_zero_pad, left=0.0, tv_core.py:326). Zero-padding makes
%   sensory_delay do double duty: it delays evidence AND deletes the first L
%   seconds of it.
%
%   Set p.n_pre = number of PRE-ONSET frames prepended to `signal` to feed the
%   REAL stream into that window instead. Then:
%     - the freeze-relative time axis is ((0:T-1) - n_pre)*dt, so index n_pre+1 is
%       t = 0 and the first n_pre samples carry negative times;
%     - the delay's interpolation pulls real samples for tau - L < 0, and the zero
%       fill applies ONLY beyond the supplied margin (the fallback);
%     - the leading n_pre samples are TRIMMED before returning, so `drift` is
%       always aligned to freeze onset (t = 0 at index 1) exactly as before.
%   n_pre = 0 (default) reproduces the bayes_fpe zero-pad byte-for-byte, so the
%   verified equivalence holds as the limiting case.
%
%   NB this is a DELIBERATE departure from bayes_fpe, which has no such mode.
%   With n_pre > 0 the causal filter also warms up on real pre-onset data rather
%   than being initialised at signal(1) -- another (documented) divergence.
%
%   full : optional 2nd output, for plotting only. Struct with the UNTRIMMED
%          transformed drift (.drift, [nTrials x T]) and its signed freeze-relative
%          time axis (.t, 1 x T, negative before onset), plus .n_pre.

    % ---- normalise input to rows [nTrials x T]; remember shape to restore ----
    was_col = iscolumn(signal);
    if isvector(signal)
        S = signal(:).';            % 1 x T
    else
        S = signal;                 % nTrials x T
    end
    T = size(S, 2);

    dt = getf(p, 'dt', []);

    % Freeze-relative time axis: st(n_pre+1) == 0, so the first n_pre samples are
    % pre-onset and carry NEGATIVE times. n_pre = 0 -> st = (0:T-1)*dt as before.
    n_pre = getf(p, 'n_pre', 0);
    assert(n_pre >= 0 && n_pre == round(n_pre) && n_pre < T, ...
           'signal_to_drift: p.n_pre must be an integer in [0, size(signal,2)).');
    if n_pre > 0
        assert(~isempty(dt), 'signal_to_drift: p.dt required when p.n_pre > 0.');
    end
    if ~isempty(dt)
        st = ((0:T-1) - n_pre) * dt;      % 1 x T, freeze-relative seconds
    else
        st = [];   % only the time-based steps need it, and they assert dt is set
    end

    % ---- 1. sensory delay (fractional-frame time shift) ---------------------
    lag = getf(p, 'sensory_delay', 0);
    if any(lag ~= 0)
        assert(~isempty(dt), 'signal_to_drift: p.dt required for sensory_delay.');
        mode = getf(p, 'sensory_delay_mode', 'delayed_evidence');
        S = apply_sensory_delay(S, lag, st, mode);
    end

    % ---- 2. causal first-order low-pass -------------------------------------
    tau = getf(p, 'filter_tau', 0);
    if any(tau > 0)
        assert(~isempty(dt), 'signal_to_drift: p.dt required for filter_tau.');
        S = causal_lowpass(S, tau, dt);
    end

    % ---- 3. saturating clip (linear up to threshold, then flat) -------------
    clip_thr = getf(p, 'clip_threshold', Inf);
    if any(isfinite(clip_thr))
        S = min(S, max(clip_thr, 0));
    end

    % ---- 4. ReLU: max(s-thr,0), optional +offset, optional ^power -----------
    relu_thr = getf(p, 'relu_threshold', 0);
    relu_off = getf(p, 'relu_offset', []);
    relu_pow = getf(p, 'relu_power', []);
    if any(relu_thr ~= 0) || ~isempty(relu_off) || ~isempty(relu_pow)
        S = relu_transform(S, relu_thr, relu_off, relu_pow);
    end

    % ---- 5. sigmoid time kernel (downweights early bins) --------------------
    mid   = getf(p, 'time_kernel_midpoint', []);
    slope = getf(p, 'time_kernel_slope', []);
    if ~isempty(mid) && ~isempty(slope)
        assert(~isempty(dt), 'signal_to_drift: p.dt required for time kernel.');
        % st is freeze-relative, so the midpoint stays anchored to onset even
        % when pre-onset samples are supplied.
        kernel = 1 ./ (1 + exp(-(st - mid) .* max(slope, 0)));
        S = S .* kernel;                             % broadcasts over rows
    end

    % ---- 6. drift = bias + gain .* processed signal -------------------------
    gain = getf(p, 'gain', 1);
    bias = getf(p, 'bias', 0);
    D = bias + gain .* S;

    if nargout > 1
        full = struct('drift', D, 't', st, 'n_pre', n_pre);
    end

    % ---- trim the pre-onset margin: drift is always aligned to t = 0 --------
    D = D(:, n_pre+1 : end);

    % ---- restore input orientation ------------------------------------------
    if isvector(signal)
        drift = D(:);
        if ~was_col, drift = drift.'; end
    else
        drift = D;
    end
end

% =========================================================================
function S = apply_sensory_delay(S, lag, st, mode)
% Fractional-frame shift via linear interpolation, per row.
%   'delayed_evidence' (default) mirrors delay_tv_signals_with_zero_pad:
%       sample at t-lag, left-edge -> 0, right-edge -> last value.
%       positive lag delays the signal and zero-pads the pre-onset interval.
%   'legacy' mirrors shift_tv_signals_by_sensory_lag:
%       sample at t+lag, left-edge -> first value, right-edge -> last value.
%
%   st is the SIGNED freeze-relative time axis. When pre-onset samples are
%   supplied (n_pre > 0) st starts negative, so queries at t-lag < 0 land on real
%   data and the left-edge fill only kicks in beyond the supplied margin.
    lag = broadcast_rows(lag, size(S, 1));
    for r = 1:size(S, 1)
        s = S(r, :);
        switch lower(string(mode))
            case {"legacy", "legacy_skip_early_signal"}
                tq = st + lag(r);
                out = interp1(st, s, tq, 'linear', NaN);
                out(tq < st(1)) = s(1);
            otherwise    % delayed_evidence / start_on_freeze_with_delayed_evidence
                tq = st - lag(r);
                out = interp1(st, s, tq, 'linear', NaN);
                out(tq < st(1)) = 0;
        end
        out(tq > st(end)) = s(end);
        S(r, :) = out;
    end
end

% =========================================================================
function S = causal_lowpass(S, tau, dt)
% First-order causal low-pass: y(1)=s(1); y(k)=a*y(k-1)+(1-a)*s(k),
% a = exp(-dt/tau)  (a=0 when tau<=0). Matches filter_tv_signals_causal.
    tau = broadcast_rows(tau, size(S, 1));
    for r = 1:size(S, 1)
        a = 0;
        if tau(r) > 0, a = exp(-dt / max(tau(r), 1e-12)); end
        s = S(r, :);
        y = s;
        for k = 2:numel(s)
            y(k) = a * y(k-1) + (1 - a) * s(k);
        end
        S(r, :) = y;
    end
end

% =========================================================================
function S = relu_transform(S, thr, offset, power)
% Exact port of relu_tv_signals_threshold (tv_core.py:480).
    thr = max(thr, 0);
    relu = max(S - thr, 0);
    if isempty(offset)
        shifted = relu;
    else
        offset = max(offset, 0);
        shifted = (relu > 0) .* (relu + offset);
    end
    if isempty(power)
        S = shifted;
        return
    end
    power = max(power, 1e-3);
    pos = max(shifted, 0);
    eps_p = 1e-6;
    powered = (pos + eps_p).^power - eps_p.^power;
    S = (pos > 0) .* powered;
end

% =========================================================================
function v = getf(p, name, default)
    if isstruct(p) && isfield(p, name) && ~isempty(p.(name))
        v = p.(name);
    else
        v = default;
    end
end

% =========================================================================
function v = broadcast_rows(v, nrows)
% Accept a scalar (same value for every row) or a per-row vector.
    if isscalar(v)
        v = repmat(v, nrows, 1);
    else
        v = v(:);
        assert(numel(v) == nrows, 'signal_to_drift: per-row param length mismatch.');
    end
end
