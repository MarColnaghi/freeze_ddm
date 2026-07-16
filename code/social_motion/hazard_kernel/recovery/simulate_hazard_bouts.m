function [S, yev, tinf, mov, slo, gid, info] = simulate_hazard_bouts(bl_t, motion_cache, beta_raw, lag_fr, fps, sp)
% SIMULATE_HAZARD_BOUTS  Simulate synthetic un-freeze events under a KNOWN hazard
% kernel β(τ), driven by the REAL social-motion each fly saw, so regressor
% autocorrelation and duration/censoring structure match the real analysis. The
% output design is exactly what fit_kernels_sparse.m consumes, so recovery runs
% the same estimation code used on real data.
%
%   bl_t         : real bout table (needs .fly, .onsets, .dur_frames; .moving_flies
%                  and .sloom used as covariates if present). dur_frames is the
%                  administrative follow-up cap for each bout.
%   motion_cache : per-fly sm series (containers.Map / indexable as in the pipeline)
%   beta_raw     : (nLag × 1) ground-truth kernel over lag_fr BEFORE amplitude
%                  calibration (SHAPE only; acausal entries should be 0). Acts on
%                  z-scored sm:  η = offset + Σ_τ β(τ)·z(sm(t−τ)).
%   lag_fr       : (nLag × 1) lag grid in frames (matches the estimator's grid)
%   fps          : frames/s
%   sp           : struct with fields
%                    dt_frames, entry_fr, n_aug, seed,
%                    target_eta_sd  (rescale β so std(kernel drive)=this),
%                    target_hazard  (sets the baseline offset = logit(target)),
%                    sm_mu, sm_sd   (global z-scoring of sm)
%
% Amplitude auto-calibration: β = β_raw · target_eta_sd / std(S·β_raw) over all
% at-risk intervals, so the effect size is realistic regardless of kernel shape.
% The SHAPE of β_raw (what recovery must return) is preserved; info.beta_true is
% the calibrated ground-truth kernel to compare recovered kernels against.

    rng(sp.seed);
    nLag    = numel(lag_fr);
    back_fr = max(lag_fr);
    fwd_fr  = -min(lag_fr);
    nB      = height(bl_t);
    mov_c   = getcol(bl_t, 'moving_flies', nB);
    slo_c   = getcol(bl_t, 'sloom',        nB);
    mask_preonset = isfield(sp,'mask_preonset') && sp.mask_preonset;  % zero pre-onset lags
    if isfield(sp,'grid_anchor'), grid_anchor = sp.grid_anchor; else, grid_anchor = 'entry'; end

    % ── Pass A: per-bout interval sm-history (computed once, reused across reps) ─
    Sb = cell(nB,1);  Tb = cell(nB,1);
    for b = 1:nB
        sm  = motion_cache(bl_t.fly(b)); sm = sm(:);
        L   = numel(sm);  on = bl_t.onsets(b);  dur = bl_t.dur_frames(b);
        if dur < sp.entry_fr, continue; end
        switch grid_anchor
            case 'entry',  base = sp.entry_fr;
            case 'offset', base = dur;
            otherwise, error('grid_anchor must be ''entry'' or ''offset''.');
        end
        kk   = ceil((1-base)/sp.dt_frames) : floor((dur-base)/sp.dt_frames);
        allg = base + kk(:)*sp.dt_frames;
        grid = allg(allg >= sp.entry_fr);
        if strcmp(grid_anchor,'entry')
            if isempty(grid) || grid(end) < dur,       grid = [grid; dur];        end %#ok<AGROW>
        else
            if isempty(grid) || grid(1) > sp.entry_fr, grid = [sp.entry_fr; grid]; end %#ok<AGROW>
        end
        Srows = nan(numel(grid), nLag);  Trows = nan(numel(grid),1);  keep = false(numel(grid),1);
        for gi = 1:numel(grid)
            f = grid(gi);  t_abs = on + f - 1;
            if (t_abs - back_fr) < 1 || (t_abs + fwd_fr) > L, continue; end
            s = (sm(t_abs - lag_fr) - sp.sm_mu) / sp.sm_sd;
            if any(isnan(s)), continue; end
            if mask_preonset, s(lag_fr > f-1) = 0; end   % drop pre-onset drive
            Srows(gi,:) = s';  Trows(gi) = f / fps;  keep(gi) = true;
        end
        Sb{b} = Srows(keep,:);  Tb{b} = Trows(keep);
    end

    % ── amplitude calibration: rescale β so kernel-drive SD = target ───────────
    Sall    = cell2mat(Sb(~cellfun(@isempty, Sb)));   % drop bouts with no valid intervals
    sd_raw  = std(Sall * beta_raw(:), 'omitnan');
    scale   = sp.target_eta_sd / max(sd_raw, eps);
    beta_true = beta_raw(:) * scale;
    offset  = log(sp.target_hazard / (1 - sp.target_hazard));

    % ── Pass B: simulate events (augmentation), stop each bout at first un-freeze ─
    nrows_b = cellfun(@(x) size(x,1), Sb);
    maxrows = sp.n_aug * sum(nrows_b);
    S    = nan(maxrows, nLag);  yev = zeros(maxrows,1);  tinf = nan(maxrows,1);
    mov  = nan(maxrows,1);      slo = nan(maxrows,1);    gid  = nan(maxrows,1);
    rr = 0;  n_cens = 0;  n_tot = 0;
    for b = 1:nB
        if isempty(Sb{b}), continue; end
        eta_b = offset + Sb{b} * beta_true;                  % reused across reps
        h_b   = 1 ./ (1 + exp(-min(max(eta_b,-30),30)));
        nib   = numel(h_b);
        for rep = 1:sp.n_aug
            ev = find(rand(nib,1) < h_b, 1, 'first');         % first un-freeze
            if isempty(ev)
                last = nib;  y = zeros(nib,1);  n_cens = n_cens + 1;   % censored at cap
            else
                last = ev;   y = zeros(ev,1);   y(ev) = 1;             % un-freeze event
            end
            n_tot = n_tot + 1;
            m = last;
            S(rr+1:rr+m, :) = Sb{b}(1:m, :);
            yev(rr+1:rr+m)  = y;
            tinf(rr+1:rr+m) = Tb{b}(1:m);
            mov(rr+1:rr+m)  = mov_c(b);
            slo(rr+1:rr+m)  = slo_c(b);
            gid(rr+1:rr+m)  = b;
            rr = rr + m;
        end
    end
    S=S(1:rr,:); yev=yev(1:rr); tinf=tinf(1:rr); mov=mov(1:rr); slo=slo(1:rr); gid=gid(1:rr);

    info = struct('beta_true',beta_true, 'offset',offset, 'scale',scale, ...
        'n_events',sum(yev), 'n_intervals',rr, 'cens_frac',n_cens/max(n_tot,1), ...
        'mean_hazard',mean(yev));
end

% ════════════════════════════════════════════════════════════════════════
function v = getcol(tbl, name, nB)
    if ismember(name, tbl.Properties.VariableNames)
        v = double(tbl.(name));  v = v(:);
    else
        v = ones(nB, 1);
    end
end
