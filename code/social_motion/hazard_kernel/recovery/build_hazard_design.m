function [S, yev, tinf, mov, slo, gid] = build_hazard_design(bl, motion_cache, cfg)
% BUILD_HAZARD_DESIGN  Turn a bout table (real or synthetic) + per-fly social
% motion into the person-period design consumed by fit_kernels_sparse. Un-freeze
% EVENTS are read from the table: one event on the final interval of each
% UNCENSORED bout (the real-pipeline convention). Mirrors the design-build loop
% of hazard_kernel_sparse.m so recovery runs identical bookkeeping to real data.
%
%   bl  : table with fields fly, onsets, dur_frames, is_censored, and (optional)
%         bout_id, moving_flies, sloom.
%   cfg : struct with fps, dt_frames, entry_fr, lag_fr, sm_mu, sm_sd.

    lag_fr  = cfg.lag_fr(:);  nLag = numel(lag_fr);
    back_fr = max(lag_fr);    fwd_fr = -min(lag_fr);
    mask_preonset = isfield(cfg,'mask_preonset') && cfg.mask_preonset;  % zero pre-onset lags
    if isfield(cfg,'grid_anchor'), grid_anchor = cfg.grid_anchor; else, grid_anchor = 'entry'; end
    nB      = height(bl);
    mov_c   = getcol(bl, 'moving_flies', nB);
    slo_c   = getcol(bl, 'sloom',        nB);
    if ismember('bout_id', bl.Properties.VariableNames)
        bid = double(bl.bout_id);  bid = bid(:);
    else
        bid = (1:nB)';
    end

    max_rows = sum(ceil(bl.dur_frames / cfg.dt_frames)) + nB;
    S = nan(max_rows, nLag);  yev = zeros(max_rows,1);  tinf = nan(max_rows,1);
    mov = nan(max_rows,1);    slo = nan(max_rows,1);    gid  = nan(max_rows,1);

    r = 0;
    for b = 1:nB
        sm  = motion_cache(bl.fly(b)); sm = sm(:);  L = numel(sm);
        on  = bl.onsets(b);  dur = bl.dur_frames(b);
        if dur < cfg.entry_fr, continue; end
        switch grid_anchor
            case 'entry',  base = cfg.entry_fr;
            case 'offset', base = dur;
            otherwise, error('grid_anchor must be ''entry'' or ''offset''.');
        end
        kk   = ceil((1-base)/cfg.dt_frames) : floor((dur-base)/cfg.dt_frames);
        allg = base + kk(:)*cfg.dt_frames;
        grid = allg(allg >= cfg.entry_fr);
        if strcmp(grid_anchor,'entry')
            if isempty(grid) || grid(end) < dur,        grid = [grid; dur];         end %#ok<AGROW>
        else
            if isempty(grid) || grid(1) > cfg.entry_fr, grid = [cfg.entry_fr; grid]; end %#ok<AGROW>
        end
        for gi = 1:numel(grid)
            f = grid(gi);  t_abs = on + f - 1;
            if (t_abs - back_fr) < 1 || (t_abs + fwd_fr) > L, continue; end
            s = (sm(t_abs - lag_fr) - cfg.sm_mu) / cfg.sm_sd;
            if any(isnan(s)), continue; end
            if mask_preonset, s(lag_fr > f-1) = 0; end   % drop pre-onset drive
            r = r + 1;
            S(r,:)  = s';  tinf(r) = f/cfg.fps;  gid(r) = bid(b);
            mov(r)  = mov_c(b);  slo(r) = slo_c(b);
            if gi == numel(grid) && ~bl.is_censored(b), yev(r) = 1; end
        end
    end
    S=S(1:r,:); yev=yev(1:r); tinf=tinf(1:r); mov=mov(1:r); slo=slo(1:r); gid=gid(1:r);
end

% ════════════════════════════════════════════════════════════════════════
function v = getcol(tbl, name, nB)
    if ismember(name, tbl.Properties.VariableNames)
        v = double(tbl.(name));  v = v(:);
    else
        v = ones(nB, 1);
    end
end
