function [bl, info] = generate_accumulator_dataset(bl_t, sm_in, P)
% GENERATE_ACCUMULATOR_DATASET  Synthetic freeze bouts from a leaky accumulator
% that integrates the REAL social-motion each focal fly experienced.
%
%   [bl, info] = generate_accumulator_dataset(bl_t, sm_in, P)
%
%   For each real freeze bout b, the per-frame drift RATE is
%       drive(t) = P.mu0 + P.beta * sm_in{b}(t)
%   i.e. a baseline urge to escape (mu0) modulated by the social motion the fly
%   saw from freeze onset (sm_in{b}, scaled, kept >= 0). A leaky accumulator
%   (sim_leaky_accumulator, leak rate P.lambda = 1/tau) integrates this drive;
%   its first-passage time across P.theta is the SYNTHETIC freeze duration.
%
%   We simulate P.n_aug independent trajectories per real bout (data
%   augmentation), each keeping that bout's fly / onset / conditions, so the
%   social-motion HISTORY behind every synthetic un-freeze is the real one and
%   the recovered hazard kernel reflects the KNOWN generative leak.
%
%   bl_t  : real template table (needs .fly, .onsets; .moving_flies, .sloom
%           used as controls if present).
%   sm_in : {nB x 1} cell, per-bout per-frame social-motion drive (onset
%           forward, capped to P.maxT_fr, already scaled; NaN -> 0). Its length
%           per bout is the integration cap (a no-crossing trajectory is
%           right-censored at that length).
%   P     : struct with fields dt, mu0, beta, sigma, lambda, theta, n_aug,
%           maxT_fr, trunc_point, seed.
%
%   Returns:
%     bl   : synthetic bout table ready for RUN_HAZARD_KERNEL — one row per
%            (bout x augmentation), with simulated .dur_frames / .is_censored,
%            filtered to .dur_frames > P.trunc_point (left truncation).
%     info : .sim_durs (uncensored, truncated, < maxT), .median_fr, .cens_frac,
%            .filt_frac (fraction lost to the trunc filter), .n_bouts.
%
%   NOTE: sim_leaky_accumulator must be on the path (code/sims/simulators).
%   Leak mapping: P.lambda = 0 -> perfect accumulator (flat kernel);
%                 P.lambda = 1/tau -> leaky (exponential kernel, timescale tau);
%                 P.lambda = 1/dt  -> extrema detector (delta kernel).

    req = {'dt','mu0','beta','sigma','lambda','theta','n_aug','maxT_fr','trunc_point','seed'};
    for i = 1:numel(req)
        assert(isfield(P, req{i}), 'generate_accumulator_dataset:missingParam', ...
               'P.%s is required', req{i});
    end

    nB    = numel(sm_in);
    n_aug = P.n_aug;
    cap_b = cellfun(@numel, sm_in);          % per-bout integration cap (frames)

    DUR  = nan(nB*n_aug, 1);
    CENS = false(nB*n_aug, 1);

    % outer loop = augmentation rep, inner = bout, so rows line up 1:1 with
    % repmat(bl_t.col, n_aug, 1) below.
    row = 0;
    for rep = 1:n_aug
        for b = 1:nB
            drive = P.mu0 + P.beta * sm_in{b}(:);          % per-frame drift RATE
            seed  = P.seed + (rep-1)*1000003 + b;           % unique per (bout,rep)
            k = sim_leaky_accumulator(drive, P.theta, P.sigma, P.lambda, P.dt, seed);
            row = row + 1;
            if isnan(k)
                DUR(row)  = cap_b(b);        % never crossed -> right-censored at cap
                CENS(row) = true;
            else
                DUR(row)  = k;               % first-passage frame (1-based)
                CENS(row) = false;
            end
        end
    end

    % ── assemble the bout table consumed by run_hazard_kernel ──────────────
    mov = getcol(bl_t, 'moving_flies', nB);
    slo = getcol(bl_t, 'sloom',        nB);

    bl = table();
    bl.fly          = repmat(bl_t.fly,    n_aug, 1);
    bl.onsets       = repmat(bl_t.onsets, n_aug, 1);
    bl.moving_flies = repmat(mov,         n_aug, 1);
    bl.sloom        = repmat(slo,         n_aug, 1);
    bl.dur_frames   = DUR;
    bl.is_censored  = CENS;
    % bout_id = the REAL bout index, SHARED across this bout's augmented copies,
    % so grouped CV in run_hazard_kernel keeps replicates with identical sm
    % history in the same fold (no augmentation leakage into the held-out set).
    bl.bout_id      = repmat((1:nB)', n_aug, 1);

    % min-duration / left-truncation filter (mirror the real pipeline)
    bl = bl(bl.dur_frames > P.trunc_point, :);

    % ── diagnostics ────────────────────────────────────────────────────────
    sd = DUR(~CENS);
    sd = sd(sd > P.trunc_point & sd < P.maxT_fr);
    info.sim_durs  = sd;
    info.median_fr = median(sd);                  % uncensored, truncated median
    info.cens_frac = mean(CENS);
    info.filt_frac = mean(DUR <= P.trunc_point);
    info.n_bouts   = height(bl);
end

% ════════════════════════════════════════════════════════════════════════
function v = getcol(tbl, name, nB)
% Return a numeric column from tbl, or ones(nB,1) if the column is absent
% (run_hazard_kernel only uses moving_flies/sloom as controls when they vary).
    if ismember(name, tbl.Properties.VariableNames)
        v = double(tbl.(name));
        v = v(:);
    else
        v = ones(nB, 1);
    end
end
